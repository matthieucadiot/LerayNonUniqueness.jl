# All cell integrals, matrix sums and checkpoints use directed Float64 intervals.
function checkpoint_intervals(A)
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),A) || error("invalid checkpoint bound")
    B=Array{Interval{Float64}}(undef,size(A))
    for i in eachindex(A)
        B[i]=interval(Float64,Float64(inf(A[i]),RoundDown),Float64(sup(A[i]),RoundUp))
    end
    return B
end

function check_matrix_column(r,row,i,nv,K)
    r["row"]==row.id && r["column"]==i && r["rows"]==collect(1:nv) || error("wrong matrix column coverage")
    size(r["G"])==size(r["G2"])==(K,K) || error("wrong matrix column dimensions")
    length(r["B2"])==row.n || error("wrong scalar norm column dimensions")
    all(A->all(x->x isa Interval{Float64}&&isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),A),(r["G"],r["G2"],r["B2"])) || error("invalid saved interval")
    return true
end

function assemble_matrix_row(inputs,B,row)
    U,V=matrix_grids(inputs,B,row)
    modes=gram_modes(row);K=length(modes)
    G=fill(I"0",K,K);G2=copy(G)
    B2=fill(I"0",row.n);columns=Vector{Int}[]
    for i in 1:length(U)-1
        file=joinpath(OUTPUT,"matrix_columns",row.id,"column_$(i).jld2")
        r=read_proof(file)
        check_matrix_column(r,row,i,length(V)-1,K)
        G=G+interval.(Float64,r["G"]);G2=G2+interval.(Float64,r["G2"])
        row.s==1 && (B2=B2+interval.(Float64,r["B2"]))
        push!(columns,[r["column"]])
    end
    check_column_coverage(columns,length(U)-1)
    if row.s==1
        B2=B2+[scalar_B2_tail(f,inputs.coercivity.bases) for f in row.fs]
    else
        B2=vector_B2(row)
    end
    br=inputs.coercivity.bases[row.r_basis];bz=inputs.coercivity.bases[row.z_basis]
    Λ=[gram_eigenvalue(br,m)+gram_eigenvalue(bz,k) for (m,k,a) in modes]
    N=assemble_N(row.d,Λ,G,G2,row.Λstar,row.W)
    design=complete_coercivity_design(row,N)
    certificate=certify_matrix(N,design.X,design.L,design.nu,row.η)
    save_proof(joinpath(OUTPUT,"matrices",row.id*".jld2");row=row,G=G,G2=G2,Λ=Λ,N=N,
               design=design,certificate=certificate,B2=B2,U=U,V=V,cells=(length(U)-1)*(length(V)-1))
    println("(v),(vi) ",row.id," verified")
    return true
end

function run_matrix_columns(inputs,B,row,first_column,last_column;deadline=Inf)
    U,V=matrix_grids(inputs,B,row)
    1≤first_column≤last_column<length(U) || error("invalid matrix column range")
    K=length(gram_modes(row));pending=Int[]
    for i in first_column:last_column
        file=joinpath(OUTPUT,"matrix_columns",row.id,"column_$(i).jld2")
        if isfile(file)
            check_matrix_column(read_proof(file),row,i,length(V)-1,K)
        else
            push!(pending,i)
        end
    end
    function completed(C)
        C.cells==length(V)-1 || error("incomplete matrix integration")
        file=joinpath(OUTPUT,"matrix_columns",row.id,"column_$(C.column).jld2")
        save_proof(file;row=row.id,column=C.column,rows=collect(1:length(V)-1),G=checkpoint_intervals(C.G),
                   G2=checkpoint_intervals(C.G2),B2=checkpoint_intervals(C.B2))
        println("saved ",row.id," column ",C.column," / ",length(U)-1);flush(stdout)
    end
    # A batch has at most one column per thread. Its vertical moments are shared;
    # each completed column is saved before starting another batch.
    for batch in Iterators.partition(pending,max(1,Threads.nthreads()))
        time()<deadline || return false
        gram_columns(inputs.profile,row,inputs.coercivity.bases,B,inputs.far,U,V,batch;on_column=completed)
    end
    return true
end
