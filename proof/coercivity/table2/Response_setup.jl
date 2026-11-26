include(joinpath(@__DIR__,"..","Paths.jl"))
# Common setup for the response scripts. Run with the Project.toml in this folder.
using IntervalArithmetic,RadiiPolynomial,LinearAlgebra,JLD2,SHA
include(joinpath(@__DIR__,"..","Compact_data.jl"))
IntervalArithmetic.configure(matmul=:slow)
BLAS.set_num_threads(1)

# Evaluate the original directed conversions once, before the response columns
# are threaded. Keep the Float64 interval and the unit-interval product separate
# so each call site retains its original interval decoration and guarantee.
const RESPONSE_PI64, RESPONSE_UNIT_PI64 = setprecision(BigFloat, 128) do
    interval(Float64, π), I"1" * π
end
const RESPONSE_BESSEL_PI, RESPONSE_BESSEL_EULER = setprecision(BigFloat, 1024) do
    interval(BigFloat, π), interval(BigFloat, Base.MathConstants.eulergamma)
end
# IntervalArithmetic adds 32 bits when reducing 1024-bit axial sin/cos arguments.
# Initialize the two directed MPFR values before those reductions run concurrently.
const RESPONSE_TRIG_PI_BOUNDS = (BigFloat(π, RoundDown; precision=1056),
                                 BigFloat(π, RoundUp; precision=1056))

include("clean_inputs.jl")
for file in ("Taylor.jl","Bsplines.jl","Residuals.jl","Asymptotics.jl","Norms.jl","Data.jl",
             "Functionals.jl","Response_data.jl","Response_fields.jl","Response_moments.jl","Response_norms.jl",
             "Response_assemble.jl","Profile_response_bounds.jl","Response_integrals.jl")
    include(file)
end

function response_parameters()
    return (prec=53,input_prec=128,N_in=100,N_band=10,N_response=16,N_ρ=32,N_functional=24,dmom=200,R0=900)
end

function table2_signature(folder,datafolder)
    files=sort(filter(f->endswith(f,".jl"),readdir(folder;join=true)))
    append!(files,[joinpath(folder,"..",f) for f in ("Project.toml","Manifest.toml","Paths.jl","Certificates.jl","Compact_data.jl")])
    result=Dict(relpath(f,dirname(folder))=>file_hash(f) for f in files)
    for name in ("responses_odd.jld2","responses_even.jld2","functionals_odd.jld2","functionals_even.jld2")
        result[name]=file_hash(joinpath(datafolder,name))
    end
    result["Profile/Interval_data.jld2"]=file_hash(joinpath(datafolder,"Profile","Interval_data.jld2"))
    compact=joinpath(COMPACT_DATA,"Profile_float64.jld2")
    isfile(compact) && (result["compact/Profile_float64.jld2"]=file_hash(compact))
    return result
end

function save_response_column(file,S,i,nv,par,signature;kind,index=0)
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&x isa Interval{Float64} && precision(inf(x))==precision(sup(x))==par.prec,S) || error("invalid column quantities")
    mkpath(dirname(file));temporary,io=mktemp(dirname(file));close(io)
    try
        jldsave(temporary;S=S,columns=[i],rows=collect(1:nv),par=par,signature=signature,kind=kind,index=index)
        mv(temporary,file)
    finally
        isfile(temporary)&&rm(temporary)
    end
end

function read_response_column(file,i,nv,par,signature;kind,index=0,nquantities=8)
    r=load(file)
    r["columns"]==[i] && r["rows"]==collect(1:nv) && r["par"]==par && r["signature"]==signature && r["kind"]==kind && r["index"]==index || error("$file does not match the requested full column")
    size(r["S"])==(nquantities,) || error("$file has incorrect dimensions")
    for (k,x) in enumerate(r["S"])
        isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&x isa Interval{Float64} && precision(inf(x))==precision(sup(x))==par.prec || error("invalid saved quantity")
        (kind=="even"&&k==6) || inf(x)≥0 || error("a squared norm is negative")
    end
    return r["S"]
end
