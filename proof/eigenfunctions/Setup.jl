using RadiiPolynomial, IntervalArithmetic, Arblib, LinearAlgebra, JLD2, SHA
IntervalArithmetic.configure(matmul=:slow)
BLAS.set_num_threads(1)
setprecision(128)
include("Functions.jl")

const DATA=get(ENV,"PROOF_DATA",joinpath(@__DIR__,"..","data"))
const OUTPUT=get(ENV,"PROOF_OUTPUT",joinpath(@__DIR__,"results"))
const N=200
const ν=interval(BigFloat,6)/interval(BigFloat,5)
const X=ℓ¹(GeometricWeight(ν))
const rstar=parse(Interval{BigFloat},"1e-20")
const WEIGHTS=[("T0",["0.5","1","15","17"],true,19), ("T1",["15","17","29","31"],true,21),
               ("T2",["29","31","42","45"],true,6), ("S0",["0.1","0.2","15","17"],true,18),
               ("S2",["29","31","152","160"],true,9), ("U0",["-7","-5","5","7"],false,18),
               ("U1",["5","7","20","23"],false,8), ("V0",["-40","-35","35","40"],false,84),
               ("V1",["-200","-180","180","200"],false,24)]
const MODES=[(name,m) for (name,nodes,radial,M) in WEIGHTS for m in 1:M]
mode_file(name,m)=name*"_m"*lpad(m,2,'0')*".jld2"

function outward128(x::Interval{BigFloat})
    isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)) || error("invalid interval enclosure")
    return interval(BigFloat,BigFloat(inf(x),RoundDown;precision=128),BigFloat(sup(x),RoundUp;precision=128))
end

function load_normalized()
    r=load(joinpath(DATA,"Eigenfunctions_normalized.jld2"))
    r["prec"]==128 || error("the normalized data need 128-bit endpoints")
    for (name,hash) in r["functional_hashes"]
        bytes2hex(sha256(read(joinpath(DATA,name))))==hash || error("a functional input changed: $name")
    end
    return r
end

function proof_signature()
    signature=Dict(name=>bytes2hex(sha256(read(joinpath(@__DIR__,name)))) for name in
                   ("Project.toml","Manifest.toml","Functions.jl","Setup.jl","Proof.jl","Run_modes.jl","Assemble.jl"))
    signature["Eigenfunctions_normalized.jld2"]=bytes2hex(sha256(read(joinpath(DATA,"Eigenfunctions_normalized.jld2"))))
    for (name,m) in MODES
        file=joinpath("Eigenfunctions",mode_file(name,m))
        signature[file]=bytes2hex(sha256(read(joinpath(DATA,file))))
    end
    return signature
end

function atomic_save(file;values...)
    mkpath(dirname(file));temporary,io=mktemp(dirname(file));close(io)
    try
        jldsave(temporary;run_id=ENV["CLEAN_RUN_ID"],values...)
        Base.Filesystem.rename(temporary,file)
    finally
        isfile(temporary)&&rm(temporary)
    end
end

function check_mode(c,name,m;fixed=nothing)
    c.name==name && c.m==m && c.N==200 && c.prec==128 && c.zeros==m || error("incorrect eigenpair coverage")
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&precision(inf(x))==precision(sup(x))==128,[c.norm2,c.r,c.Z,c.epsilon,c.theta...,c.errors...]) || error("a certificate contains invalid intervals")
    inf(c.norm2)>0 && inf(c.theta[2])>0 || error("normalization or positive phase failed")
    sup(c.Z)<1 && sup(c.r)<=inf(rstar) || error("matching contraction failed")
    sup(c.r)<=inf(parse(Interval{BigFloat},"3.3e-31")) || error("the eigenvalue target failed")
    all(x->isguaranteed(x)&&isfinite(sup(x))&&inf(x)>=0&&sup(x)<=inf(parse(Interval{BigFloat},"1.1e-24")),c.errors) || error("a piecewise derivative target failed")
    isguaranteed(c.lambda) || error("unguaranteed eigenvalue")
    if fixed!==nothing
        for key in (:lambda_bar,:A,:B)
            isequal_interval(getproperty(c,key),getproperty(fixed,key)) || error("the certified target changed")
        end
        for key in (:phi_a,:psi_a,:psi_b,:phi_b)
            all(isequal_interval.(getproperty(c,key),getproperty(fixed,key))) || error("the certified polynomial changed")
        end
    end
    return true
end
