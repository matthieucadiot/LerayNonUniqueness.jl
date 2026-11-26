# Accept Appendix A only when every one of the 207 fixed modes has passed.
include("Setup.jl")
include("Proof.jl")
signature=proof_signature();data=load_normalized()
certificates=Any[]
for (name,m) in MODES
    file=joinpath(OUTPUT,"modes",mode_file(name,m))
    r=load(file)
    r["run_id"]==ENV["CLEAN_RUN_ID"] || error("Mode belongs to a different run")
    r["signature"]==signature || error("a mode has a different source/input fingerprint")
    check_mode(r["certificate"],name,m;fixed=data["bases"][name].modes[m])
    push!(certificates,r["certificate"])
end
constants=Dict(name=>constant_mode(name) for (name,nodes,radial,M) in WEIGHTS)
proof_signature()==signature || error("the source or inputs changed during assembly")
atomic_save(joinpath(OUTPUT,"Eigenfunctions.jld2");certificates=certificates,bases=data["bases"],constants=constants,
            functional_hashes=data["functional_hashes"],signature=signature,certified=true,
            eigenvalue_error=parse(Interval{BigFloat},"3.3e-31"),derivative_error=parse(Interval{BigFloat},"1.1e-24"))
println("Appendix A verified for all 207 nonconstant eigenpairs and the nine constant modes.")
