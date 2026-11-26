# One shared sweep for the residual table and every bound in Table 1.
started=time_ns()
include("Setup.jl")
include("Files.jl")
length(ARGS)==2 || error("usage: Run_columns.jl FIRST LAST")
first=parse(Int,ARGS[1]);last=parse(Int,ARGS[2])
seconds=parse(Int,get(ENV,"PROOF_JOB_SECONDS","3300"));seconds≥0 || error("invalid soft time budget")
expired()=(time_ns()-started)/1_000_000_000≥seconds
par=residual_parameters();setprecision(par.prec)
data=load_data(joinpath(DATA,"Profile"))
nu=length(data.U)-1;nv=length(data.V)-1
1≤first≤last≤nu || error("invalid column range")
signature=residual_signature();pending=Int[];prefix=Dict{Int,NamedTuple}()
for i in first:last
    if isfile(column_file(i))
        read_column(column_file(i),i,nv,par,signature)
        isfile(partial_file(i))&&rm(partial_file(i))
    else
        prefix[i]=isfile(partial_file(i)) ? read_column(partial_file(i),i,nv,par,signature;complete=false) :
                                         (S=zeros(Interval{BigFloat},19),rows_done=0)
        push!(pending,i)
    end
end
if !isempty(pending)&&!expired()
    far=far_field_data(par.prec,data.λ̃)
    row=[row_data(j,data,par) for j in 1:nv]
    writing=ReentrantLock()
    println("Residual and Table 1: columns $first:$last, $nv rows, $(Threads.nthreads()) threads")
    flush(stdout)
    Threads.@threads :dynamic for i in pending
        setprecision(par.prec) do
            if !expired()
                C=column_data(i,data,par);s=prefix[i].S;done=prefix[i].rows_done
                for j in done+1:nv
                    expired()&&break
                    accumulate!(s,profile_cell(C,row[j],data,far,par));done=j
                    if done%16==0&&done<nv
                        lock(writing) do
                            save_column(partial_file(i),s,i,nv,done,par,signature)
                        end
                    end
                end
                lock(writing) do
                    if done==nv
                        save_column(column_file(i),s,i,nv,done,par,signature)
                        isfile(partial_file(i))&&rm(partial_file(i))
                        println("column $i complete")
                    else
                        save_column(partial_file(i),s,i,nv,done,par,signature)
                        println("column $i checkpointed through row $done/$nv")
                    end
                    flush(stdout)
                end
            end
        end
    end
end
residual_signature()==signature || error("the input or numerical code changed during the run")
complete=all(isfile(column_file(i)) for i in first:last)
println(complete ? "Range complete." : "Soft deadline reached; saved prefixes will resume.")
println("Peak RSS: $(round(Sys.maxrss()/2.0^30;digits=2)) GiB")
exit(complete ? 0 : 75)
