# Store only the fixed normalized interval approximations, without old proof bounds.
using IntervalArithmetic, RadiiPolynomial, JLD2, SHA
setprecision(128)
folder=get(ENV,"PROOF_DATA",joinpath(@__DIR__,"..","data"))
prop=load(joinpath(folder,"Proposition_4_21_inputs.jld2"),"bases")
bases=Dict{String,Any}()
for (name,b) in prop
    modes=[(m=m.m,lambda_bar=m.lambda_bar,phi_a=m.phi_a,psi_a=m.psi_a,
            psi_b=m.psi_b,phi_b=m.phi_b,A=m.A,B=m.B) for m in b.modes]
    bases[name]=(name=b.name,radial=b.radial,nodes=b.nodes,N=b.N,constant=b.constant,modes=modes)
end
hashes=Dict{String,String}()
for filename in ("functionals_odd.jld2","functionals_even.jld2")
    file=joinpath(folder,filename)
    hashes[filename]=bytes2hex(sha256(read(file)))
    for (name,b) in load(file,"eigenfunctions"), mode in b.modes
        fixed=bases[name].modes[mode.m]
        for key in (:lambda_bar,:A,:B)
            isequal_interval(getproperty(mode,key),getproperty(fixed,key)) || error("the functional target differs: $name $(mode.m)")
        end
        for key in (:phi_a,:psi_a,:psi_b,:phi_b)
            all(isequal_interval.(getproperty(mode,key),getproperty(fixed,key))) || error("a functional polynomial differs: $name $(mode.m)")
        end
    end
end
length(bases)==9 && sum(length(b.modes) for b in values(bases))==207 || error("the nine bases and 207 modes are required")
# This preparation is an authoring step; the proof loads this JLD2 file directly.
jldsave(joinpath(folder,"Eigenfunctions_normalized.jld2");bases=bases,functional_hashes=hashes,prec=128)
println("Saved 207 fixed normalized interval approximations, exactly matching the functional inputs.")
