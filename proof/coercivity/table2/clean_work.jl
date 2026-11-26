include("Response_setup.jl")
include("clean_geometry.jl")
using TOML
const CLEAN_COLUMN_IO = ReentrantLock()
const CLEAN_SCOPE_SHA = "0837f60b3181ba9da9da12f1025f569c24cccccd329edda9bb6d21521c7d9159"

function clean_work_plan()
    scope=joinpath(@__DIR__,"repair_scope.tsv")
    clean_hash(scope)==CLEAN_SCOPE_SHA || error("Unexpected junction-refinement column list")
    repairs=Set{Tuple{Bool,Int,Int}}()
    lines=readlines(scope)
    first(lines)=="lane\tkind\tindex\tcolumn" || error("Invalid column list")
    for row in lines[2:end]
        _,kind,index,column=split(row,'\t')
        kind in ("odd","even") || error("Invalid parity in column list")
        push!(repairs,(kind=="odd",parse(Int,index),parse(Int,column)))
    end
    length(repairs)==length(lines)-1==4289 || error("Incomplete or duplicate refinement inventory")
    return repairs
end

function clean_run_id()
    id=get(ENV,"CLEAN_RUN_ID","")
    occursin(r"^[0-9a-f]{32}$",id) || error("CLEAN_RUN_ID must be the 32-digit hexadecimal ID of this run")
    return id
end

function clean_signature()
    signature=table2_signature(@__DIR__,PROOF_DATA)
    signature["clean-run-id"]=clean_run_id()
    signature["refinement-scope"]=CLEAN_SCOPE_SHA
    signature["fixed-coefficients"]=clean_hash(joinpath(PROOF_DATA,"Coercivity",CLEAN_INPUT_NAME))
    signature["weights"]=CLEAN_WEIGHTS_SHA
    # Every table2 Julia source, all four source datasets, and the arithmetic
    # project and fixed profile enter the column identity.
    io=IOBuffer()
    for key in sort!(collect(keys(signature)))
        println(io,key,"\t",signature[key])
    end
    bytes2hex(sha256(take!(io)))
end

function clean_geometry(t,odd,index,inventory)
    data=response_profile(t,odd); qs=odd ? t.qo : t.qe; f=qs.functionals[index]
    U=refinement_axis_plan(data.U,refinement_nodes(f,qs.eigenfunctions,:u))
    V=refinement_axis_plan(data.V,refinement_nodes(f,qs.eigenfunctions,:v))
    active=if odd && index in 1:9
        [48]
    elseif !odd && index in 1:12
        [49]
    elseif (odd && index in 22:24) || (!odd && index in 22:32)
        [i for i in 1:length(data.U)-1 if begin
            I=refinement_physical_interval(data.U[i:i+1])
            sup(I)>inf(f.r_spec[1]) && inf(I)<sup(f.r_spec[4])
        end]
    else
        Int[]
    end
    exact=Set(i for (o,j,i) in inventory if o==odd && j==index)
    exact==Set(active) || error("Geometric scope differs for response $odd/$index")
    return (U=U,V=V,columns=exact)
end

clean_column_path(odd,j,i)=joinpath(COERCIVITY_RESULTS,"table2_columns",(odd ? "odd_" : "even_")*string(j),"column_$i.jld2")
clean_partial_path(path)=replace(path,r"column_(\d+)\.jld2$"=>s"partial_\1.jld2")

function clean_save(writer,path)
    mkpath(dirname(path));tmp,io=mktemp(dirname(path));close(io)
    try
        writer(tmp)
        Base.Filesystem.rename(tmp,path)
    finally
        isfile(tmp) && rm(tmp)
    end
end

function clean_children(plan,i,j)
    i in plan.columns || return 1
    (length(plan.U.partitions[i])-1)*(length(plan.V.partitions[j])-1)
end
function clean_expected_children(plan,i,n)
    i in plan.columns || return n
    (length(plan.U.partitions[i])-1)*sum(j->length(plan.V.partitions[j])-1,1:n;init=0)
end

function clean_read(file,odd,index,col,nv,plan,signature;complete)
    r=lock(CLEAN_COLUMN_IO) do
        load(file)
    end
    r["schema"]=="clean-weighted-column-v1" && r["signature"]==signature &&
        r["odd"]===odd && r["index"]==index && r["column"]==col && r["total_rows"]==nv &&
        r["refined"]==in(col,plan.columns) || error("Wrong column identity: $file")
    n=r["rows_done"]
    n isa Int && 0<=n<=nv && r["complete"]==complete &&
        (complete ? n==nv : n<nv) || error("Incomplete or invalid column: $file")
    r["children"]==clean_expected_children(plan,col,n) || error("Child coverage mismatch: $file")
    return (S=refinement_validate_quantities(r["S"];odd), rows=n, children=r["children"])
end

function clean_write(file,S,odd,index,col,nv,plan,signature,rows,children)
    refinement_validate_quantities(S;odd)
    children==clean_expected_children(plan,col,rows) || error("Incomplete child-cell sum")
    lock(CLEAN_COLUMN_IO) do
        clean_save(file) do tmp
            jldsave(tmp;schema="clean-weighted-column-v1",signature,odd,index,column=col,
                    total_rows=nv,refined=in(col,plan.columns),rows_done=rows,children,complete=rows==nv,S)
        end
    end
end
