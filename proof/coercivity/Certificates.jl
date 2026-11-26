using IntervalArithmetic,RadiiPolynomial,JLD2,SHA
include("Paths.jl")

function check_residual_certificate(r,par,profile_hash)
    r["certified"]===true && r["columns"]==collect(1:804) && r["rows"]==collect(1:575) || error("incomplete residual certificate")
    size(r["S"])==(19,804) || error("wrong residual column dimensions")
    all(x->x isa Interval{BigFloat}&&isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&inf(x)≥0&&
           precision(inf(x))==precision(sup(x))==128,r["S"]) || error("invalid residual interval columns")
    r["profile_hash"]==profile_hash || error("the residual certificate uses different profiles")
    r["par"].prec==128 || error("the residual certificate must use 128-bit arithmetic")
    for key in (:N_in,:N_band,:N_ρ,:N_functional,:dmom,:R0)
        getproperty(r["par"],key)==getproperty(par,key) || error("different profile parameter: $key")
    end
    c=r["profile_constants"]
    c.certified===true || error("the profile constants are conditional")
    for key in (:Cb,:RΨZ,:sΩ,:IoΩ,:sΩfar,:rRΨ̃,:gΨ̃_min,:nΩ̃Z,:sSΩ̃,:ℛ̃,:δ)
        x=getproperty(c,key)
        x isa Interval{BigFloat}&&isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&inf(x)≥0 || error("invalid profile constant: $key")
    end
    isthin(c.gΨ̃_min)&&inf(c.gΨ̃_min)>0 || error("invalid lower bound for the even gradient")
    in_interval(big(2)//big(25),c.Cb)&&in_interval(big(1)//big(10)^6,c.δ) || error("the fixed constants differ")
    return c
end

function load_residual_profile_constants(par)
    file=get(ENV,"RESIDUAL_CERTIFICATE",normpath(joinpath(@__DIR__,"..","residual","results","Residual_Table1.jld2")))
    isfile(file) || error("complete the residual/Table 1 proof before final Table 2 assembly: $file")
    r=load(file)
    r["run_id"]==ENV["CLEAN_RUN_ID"] || error("Upstream certificate belongs to a different run: $file")
    c=check_residual_certificate(r,par,bytes2hex(sha256(read(joinpath(PROOF_DATA,"Profile","Interval_data.jld2")))))
    for (name,hash) in r["signature"]
        source=name=="Profile/Interval_data.jld2" ? joinpath(PROOF_DATA,name) : joinpath(@__DIR__,"..","residual",name)
        hash==bytes2hex(sha256(read(source))) || error("a residual input or source changed: $name")
    end
    # The residual is certified at 128 bits; the coercivity kernels use outward Float64 enclosures.
    c=map(x->x isa Interval ? interval(Float64,x) : x,c)
    return c,(file=basename(file),sha256=bytes2hex(sha256(read(file))),signature=r["signature"])
end

function check_shared_sources()
    p=load(joinpath(PROOF_DATA,"Proposition_4_21_inputs.jld2"),"sources")
    files=("Interval_data.jld2"=>joinpath(PROOF_DATA,"Profile","Interval_data.jld2"),
           "functionals_odd.jld2"=>joinpath(PROOF_DATA,"functionals_odd.jld2"),
           "functionals_even.jld2"=>joinpath(PROOF_DATA,"functionals_even.jld2"))
    for (name,file) in files
        p[name]==bytes2hex(sha256(read(file))) || error("Table 2 and Proposition 4.21 use different $name")
    end
    return true
end

function load_spectral_certificate()
    file=get(ENV,"EIGENFUNCTIONS_CERTIFICATE",normpath(joinpath(@__DIR__,"..","eigenfunctions","results","Eigenfunctions.jld2")))
    isfile(file) || error("complete the eigenfunction proof before final coercivity assembly: $file")
    r=load(file)
    r["run_id"]==ENV["CLEAN_RUN_ID"] || error("Upstream certificate belongs to a different run: $file")
    r["certified"]===true && length(r["certificates"])==207 || error("incomplete Appendix A certificate")
    for name in ("functionals_odd.jld2","functionals_even.jld2")
        r["functional_hashes"][name]==bytes2hex(sha256(read(joinpath(PROOF_DATA,name)))) || error("the spectral certificate uses different functionals")
    end
    bases=load(joinpath(PROOF_DATA,"Proposition_4_21_inputs.jld2"),"bases")
    constants=r["constants"]
    Set(keys(constants))==Set(keys(bases)) && length(constants)==9 || error("Incomplete constant-mode certificate")
    all(name->issubset_interval(constants[name],bases[name].constant),keys(bases)) ||
        error("A fixed constant mode does not enclose its certified value")
    expected=Set((name,m) for (name,b) in bases for m in 1:length(b.modes))
    Set((c.name,c.m) for c in r["certificates"])==expected || error("the eigenpair coverage differs")
    for c in r["certificates"]
        c.N==200&&c.prec==128&&c.zeros==c.m || error("incorrect spectral indexing")
        all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),(c.norm2,c.r,c.Z,c.lambda)) || error("invalid spectral bounds")
        inf(c.norm2)>0&&inf(c.r)≥0&&sup(c.r)≤inf(big(I"3.3e-31"))&&sup(c.Z)<1 || error("a spectral target failed")
        length(c.errors)==3&&all(x->isguaranteed(x)&&isfinite(sup(x))&&inf(x)≥0&&sup(x)≤inf(big(I"1.1e-24")),c.errors) || error("a mode error exceeds Appendix A")
        fixed=bases[c.name].modes[c.m]
        issubset_interval(c.lambda,fixed.lambda) || error("A fixed eigenvalue does not enclose its certified value")
        for key in (:lambda_bar,:A,:B)
            isequal_interval(getproperty(c,key),getproperty(fixed,key)) || error("the certified mode differs from the coercivity input")
        end
        for key in (:phi_a,:psi_a,:psi_b,:phi_b)
            all(isequal_interval.(getproperty(c,key),getproperty(fixed,key))) || error("the certified mode coefficients differ")
        end
    end
    for (name,hash) in r["signature"]
        source=endswith(name,".jld2") ? joinpath(PROOF_DATA,name) : joinpath(@__DIR__,"..","eigenfunctions",name)
        hash==bytes2hex(sha256(read(source))) || error("an eigenfunction input or source changed: $name")
    end
    return (file=basename(file),sha256=bytes2hex(sha256(read(file))),signature=r["signature"])
end
