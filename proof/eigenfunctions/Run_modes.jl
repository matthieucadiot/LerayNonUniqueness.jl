# Usage: julia --project=. Run_modes.jl FIRST LAST. Each accepted mode is atomic.
started=time_ns()
include("Setup.jl")
include("Proof.jl")
function main()
    length(ARGS)==2 || error("usage: Run_modes.jl FIRST LAST")
    first,last=parse.(Int,ARGS)
    1<=first<=last<=207 || error("the mode range must lie in 1:207")
    seconds=parse(Int,get(ENV,"PROOF_JOB_SECONDS","3300"));seconds>=0 || error("invalid soft time budget")
    signature=proof_signature();data=load_normalized();roots=roots_Shat()
    contexts=Dict{String,Any}()
    complete=true
    for i in first:last
        name,m=MODES[i];file=joinpath(OUTPUT,"modes",mode_file(name,m))
        if isfile(file)
            saved=load(file)
            saved["run_id"]==ENV["CLEAN_RUN_ID"] || error("Mode belongs to a different run")
            saved["signature"]==signature || error("the saved mode has a different input or source fingerprint")
            check_mode(saved["certificate"],name,m;fixed=data["bases"][name].modes[m])
            continue
        end
        if (time_ns()-started)/1e9>=seconds
            complete=false;break
        end
        if !haskey(contexts,name)
            contexts[name]=weight_context(name,roots)
        end
        result=prove_mode(contexts[name],m,data["bases"][name].modes[m])
        atomic_save(file;certificate=result,signature=signature)
        println(name," ",m,": r=",result.r,"; errors=",result.errors,"; zeros=",result.zeros);flush(stdout)
    end
    proof_signature()==signature || error("the source or inputs changed during the proof")
    println("Peak RSS bytes: ",Sys.maxrss())
    return complete ? 0 : 75
end
exit(main())
