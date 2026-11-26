# Exercise the corrected fixed inputs, both parities, a junction split and the
# worker's resumable prefix on actual mathematical cells. Output is temporary.
include("clean_work.jl")
using Test

function clean_smoke()
    par=response_parameters();setprecision(BigFloat,par.input_prec)
    file=clean_load_design(PROOF_DATA)
    # Use separate data structures to detect repeated normalization or mutation.
    clean=load_table2_data(PROOF_DATA;clean=file)
    old=load_table2_data(PROOF_DATA;clean=nothing)
    @test all(isequal_interval.(clean.even.μ,old.even.μ))
    for (j,w) in enumerate(clean_design(PROOF_DATA))
        @test isequal_interval(clean.qo.functionals[j].nu,
            interval(Float64,old.qo.functionals[j].nu)*interval(Float64,w.factor)^2)
    end
    inventory=clean_work_plan()
    stress=NamedTuple[]
    for (odd,index) in ((true,21),(true,22),(false,24))
        plan=clean_geometry(clean,odd,index,inventory)
        qs=odd ? clean.qo : clean.qe; data=response_profile(clean,odd); f=qs.functionals[index]
        active(axis,spec)=findall(i->begin
            I=refinement_physical_interval(axis[i:i+1]);sup(I)>inf(spec[1]) && inf(I)<sup(spec[4])
        end,1:length(axis)-1)
        us=active(data.U,f.r_spec);vs=active(data.V,f.z_spec)
        i=us[cld(length(us),2)];j=vs[cld(length(vs),2)]
        if i in plan.columns
            splitrows=filter(k->length(plan.V.partitions[k])>2,vs)
            isempty(splitrows) || (j=first(splitrows))
        end
        unit=(odd=odd,index=index,column=i)
        geometry=i in plan.columns ? plan : (U=(partitions=[data.U[k:k+1] for k in 1:length(data.U)-1],),
                                              V=(partitions=[data.V[k:k+1] for k in 1:length(data.V)-1],))
        newcell=refinement_evaluate_row(refinement_prepare_case(clean,par,unit,geometry),j)
        oldcell=refinement_evaluate_row(refinement_prepare_case(old,par,unit,geometry),j)
        factor=odd ? clean_design(PROOF_DATA)[index].factor : 1//1
        @test all(x->!isempty_interval(x),intersect_interval.(newcell.S,oldcell.S .* interval(Float64,factor)^2))
        @test newcell.children == oldcell.children == clean_children(plan,i,j)
        println("Checked ",odd ? "odd" : "even"," response ",index," cell ",i,",",j,
                " with ",newcell.children," child cell(s); all eight independently evaluated enclosures overlap.")
        push!(stress,(odd=odd,index=index,unit=unit,plan=geometry,row=j,expected=newcell.S))
        if !odd
            @test any(x->sup(x)>0,newcell.S)
        end
    end
    # Simultaneous 128-bit outer scopes, 256-bit spline conversions and
    # 1024-bit Bessel seeds must agree with the serial enclosure endpoints.
    results=Vector{Bool}(undef,4*length(stress))
    Threads.@threads :dynamic for k in eachindex(results)
        c=stress[mod1(k,length(stress))]
        results[k]=setprecision(BigFloat,128) do
            ctx=refinement_prepare_case(clean,par,c.unit,c.plan)
            cell=refinement_evaluate_row(ctx,c.row)
            probe=setprecision(BigFloat,256) do
                yield()
                precision(BigFloat)
            end
            probe==256 && precision(BigFloat)==128 && all(isequal_interval.(c.expected,cell.S))
        end
    end
    @test all(results)
    println(length(results)," repeated cell evaluations on ",Threads.nthreads()," threads matched all serial endpoints.")
    mktempdir() do tmp
        plan=(columns=Set{Int}(),); s=zeros(Interval{Float64},8)
        partial=joinpath(tmp,"partial_4.jld2")
        clean_write(partial,s,true,1,4,809,plan,"test",1,1)
        @test clean_read(partial,true,1,4,809,plan,"test";complete=false).rows==1
        @test_throws ErrorException clean_read(partial,true,1,4,809,plan,"other";complete=false)
        @test_throws ErrorException clean_read(partial,true,1,4,809,plan,"test";complete=true)
    end
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    @testset "Fixed data, physical cells, and column restart" clean_smoke()
end
