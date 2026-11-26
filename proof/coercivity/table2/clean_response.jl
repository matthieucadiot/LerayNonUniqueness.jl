# Run one disjoint lane through all response grids. The numerical identity is
# independent of worker count; a new process resumes an atomic row prefix.
const CLEAN_START = time_ns()
include("clean_work.jl")

function work_on_column(t,par,odd,index,col,plan,signature,caches,expired)
    data=response_profile(t,odd);nv=length(data.V)-1; out=clean_column_path(odd,index,col)
    partial=clean_partial_path(out)
    state=isfile(partial) ? clean_read(partial,odd,index,col,nv,plan,signature;complete=false) :
                            (S=zeros(Interval{Float64},8),rows=0,children=0)
    rd=odd ? t.odd : t.even; qs=odd ? t.qo : t.qe
    refined=col in plan.columns
    ctx=refined ? refinement_prepare_case(t,par,(odd=odd,index=index,column=col),plan) : nothing
    C=refined ? nothing : column_data(col,data,par)
    CR=refined ? nothing : response_basis_column(col,rd,data.U)
    S=state.S;rows=state.rows;children=state.children
    for j in rows+1:nv
        expired() && break
        if refined
            cell=refinement_evaluate_row(ctx,j)
            cell.coverage_checked && cell.children==clean_children(plan,col,j) || error("Missing child")
            S=S+cell.S;children+=cell.children
        else
            R=caches.R[j];P=profile_fields_on_cell(data,caches.far,C,R,par)
            Q=functional_on_cell(qs.functionals[index],qs.eigenfunctions,C,R,par)
            value=vcat(response_cell_quantities(index,rd,CR,caches.RR[j],C,R,P,caches.ext,Q.q,par),
                       functional_cell_integrals(Q,odd,C,R))
            S=S+refinement_validate_quantities(value;odd);children+=1
        end
        rows=j
        if rows%32==0 && rows<nv
            clean_write(partial,S,odd,index,col,nv,plan,signature,rows,children)
        end
    end
    clean_write(rows==nv ? out : partial,S,odd,index,col,nv,plan,signature,rows,children)
    if rows==nv
        isfile(partial) && rm(partial)
        println(odd ? "odd " : "even ",index," column ",col," complete");flush(stdout)
    end
end

function clean_response(lane,lanes,seconds)
    1<=lane<=lanes && seconds>0 || error("Invalid lane or soft limit")
    par=response_parameters();setprecision(BigFloat,par.input_prec)
    expired()=(time_ns()-CLEAN_START)/1.0e9 >= seconds
    inventory=clean_work_plan();signature=clean_signature()
    t=load_table2_data(PROOF_DATA;clean=clean_load_design(PROOF_DATA))
    for odd in (true,false)
        data=response_profile(t,odd);rd=odd ? t.odd : t.even
        nu=length(data.U)-1;nv=length(data.V)-1
        (nu,nv)==(odd ? (1157,809) : (1175,888)) || error("Response domain changed")
        local caches_base
        built=false
        for index in 1:rd.n
            # Paths only schedule work. The final assembly re-opens every column
            # and checks its fingerprint, row/child coverage and interval contents.
            todo=[i for i in lane:lanes:nu if !isfile(clean_column_path(odd,index,i))]
            isempty(todo) && continue
            expired() && return 75
            plan=clean_geometry(t,odd,index,inventory)
            if !built
                far=far_field_data(par.prec,data.λ̃)
                caches_base=(R=[row_data(j,data,par) for j in 1:nv],
                             RR=[response_basis_row(j,rd,data.V) for j in 1:nv],far=far)
                built=true
            end
            caches=merge(caches_base,(ext=response_farfield(index,rd,caches_base.far),))
            Threads.@threads :dynamic for col in todo
                if !expired()
                    setprecision(BigFloat,par.input_prec) do
                        work_on_column(t,par,odd,index,col,plan,signature,caches,expired)
                    end
                end
            end
            all(i->isfile(clean_column_path(odd,index,i)),todo) || return 75
        end
    end
    return 0
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==3 || error("usage: clean_response.jl LANE LANES SOFT_SECONDS")
    exit(clean_response(parse.(Int,ARGS)...))
end
