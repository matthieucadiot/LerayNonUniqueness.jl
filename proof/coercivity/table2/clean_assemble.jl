# Assemble Table 2 and Proposition 4.21(vii) from this run's complete columns.
include("clean_work.jl")
include(joinpath(@__DIR__,"..","Certificates.jl"))

function clean_assemble()
    par=response_parameters();setprecision(BigFloat,par.input_prec)
    signature=clean_signature();inventory=clean_work_plan()
    t=load_table2_data(PROOF_DATA;clean=clean_load_design(PROOF_DATA))
    check_shared_sources() || error("Different raw spectral data")
    eigen=load_spectral_certificate()
    c,certificate=load_residual_profile_constants(par)
    groups=Dict{Bool,Any}();physical=Dict{String,Any}()
    tails=Dict{Int,Interval{Float64}}();column_sums=Dict{Tuple{Bool,Int},Vector{Interval{Float64}}}()
    counts=Dict("columns"=>0,"selected_columns"=>0,"parent_cells"=>0,"evaluated_cells"=>0)
    radius=sqrt(U0(interval(Float64,t.original_profile.L)))
    for odd in (true,false)
        side=odd ? "odd" : "even"
        rd=odd ? t.odd : t.even; qs=odd ? t.qo : t.qe
        data=response_profile(t,odd);nu=length(data.U)-1;nv=length(data.V)-1
        (nu,nv,rd.n)==(odd ? (1157,809,26) : (1175,888,32)) || error("Invalid coverage")
        far=far_field_data(par.prec,data.λ̃)
        raw=NamedTuple[];functionals=NamedTuple[];squared=Interval{Float64}[]
        for index in 1:rd.n
            plan=clean_geometry(t,odd,index,inventory)
            S=zeros(Interval{Float64},8)
            for i in 1:nu
                record=clean_read(clean_column_path(odd,index,i),odd,index,i,nv,plan,signature;complete=true)
                S=S+record.S
                counts["columns"]+=1
                counts["selected_columns"]+=i in plan.columns
                counts["parent_cells"]+=nv
                counts["evaluated_cells"]+=record.children
            end
            refinement_validate_quantities(S;odd)
            column_sums[(odd,index)]=S
            f=qs.functionals[index];scalar=f.kind=="scalar"
            scalar==(!odd && 24<=index<=32) || error("Unexpected scalar source")
            tail=zero(S[8])
            if scalar
                rows=refinement_tail_rows(t,par,index,plan)
                tail=functional_tail_integral(f,qs.eigenfunctions,rows.refined,par,interval(Float64,sup(f.r_spec[4])))
                response_nonnegative(tail) || error("Invalid scalar tail $index")
                tails[index]=tail
            end
            exterior=odd ? nothing : exterior_response_quantities(response_farfield(index,rd,far),radius)
            push!(raw,raw_response_bounds(S,odd,exterior;complete=true))
            norm2=S[8]+tail;push!(squared,norm2)
            push!(functionals,(norm=min(sqrt(S[7]),sqrt(norm2)),error=zero(S[8]),certified=qs.published_mode_bounds))
            println("Validated ",side," response ",index," (",nu," columns).");flush(stdout)
        end
        roots=sum(sqrt.(squared));squares=sum(squared)
        rlimit,slimit=odd ? (96//5,993//50) : (34//1,61//1)
        all(response_nonnegative,squared) && sup(roots)<=inf(interval(Float64,rlimit)) &&
            sup(squares)<=inf(interval(Float64,slimit)) || error("Proposition 4.21(vii) failed for $side")
        physical[side]=Dict("sum_B_upper"=>sup(roots),"sum_B_squared_upper"=>sup(squares),
                            "individual_B_squared_upper"=>sup.(squared))
        groups[odd]=(raw=raw,functionals=functionals)
    end
    counts["columns"]==67682 && counts["selected_columns"]==4289 && counts["parent_cells"]==57725138 ||
        error("Missing columns or cells")
    Set(keys(tails))==Set(24:32) || error("Missing scalar tails")
    table=assemble_responses(groups[true].raw,groups[false].raw,groups[true].functionals,groups[false].functionals,c)
    verify_table2(table)
    result=Dict{String,Any}("all_checks_passed"=>true,"physical_passed"=>true,"schema"=>"clean-table2-v1",
        "identity"=>signature,"run_id"=>clean_run_id(),
        "fixed_coefficients_sha256"=>clean_hash(joinpath(PROOF_DATA,"Coercivity",CLEAN_INPUT_NAME)),
        "scope_sha256"=>CLEAN_SCOPE_SHA,"weights_sha256"=>CLEAN_WEIGHTS_SHA,"coverage"=>counts,"physical"=>physical,"scalar_tail_upper"=>Dict(string(j)=>sup(v) for (j,v) in tails),
        "eigenfunctions_sha256"=>eigen.sha256,"residual_sha256"=>certificate.sha256)
    for side in (:odd,:even)
        result[string(side)]=Dict(string(k)=>Dict("upper"=>sup(x),"target_lower"=>inf(table2_targets()[side][k]),
            "passed"=>true) for (k,x) in pairs(table[side]))
    end
    clean_save(joinpath(COERCIVITY_RESULTS,"Table2_bounds.jld2")) do path
        jldsave(path;table,physical,scalar_tails=tails,column_sums,signature,counts,run_id=clean_run_id(),
                profile_certificate=certificate,eigenfunctions_sha256=eigen.sha256)
    end
    result["table2_data_sha256"] = clean_hash(joinpath(COERCIVITY_RESULTS,"Table2_bounds.jld2"))
    clean_save(joinpath(COERCIVITY_RESULTS,"table2_checks.toml")) do path
        open(io->TOML.print(io,result;sorted=true),path,"w")
    end
    println("All eight Table 2 targets and both physical norm sums per parity passed.")
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    clean_assemble()
end
