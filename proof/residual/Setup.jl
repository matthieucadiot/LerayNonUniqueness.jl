# Common arithmetic and parameters for the residual and all Table 1 bounds.
using IntervalArithmetic,RadiiPolynomial,LinearAlgebra,JLD2,SHA
IntervalArithmetic.configure(matmul=:slow)
BLAS.set_num_threads(1)
for file in ("Taylor.jl","Bsplines.jl","Residuals.jl","Asymptotics.jl","Norms.jl","Data.jl",
             "Response_moments.jl","Profile_integrals.jl","Profile_response_bounds.jl","Table1_fields.jl")
    include(file)
end

const DATA=get(ENV,"PROOF_DATA",joinpath(@__DIR__,"..","data"))
const OUTPUT=get(ENV,"PROOF_OUTPUT",joinpath(@__DIR__,"results"))
residual_parameters()=(prec=128,N_in=100,N_band=10,N_ρ=32,N_functional=24,dmom=200,R0=900)
profile_file()=joinpath(DATA,"Profile","Interval_data.jld2")
file_hash(file)=open(file) do io
    bytes2hex(sha256(io))
end

function residual_signature()
    names=("Setup.jl","Taylor.jl","Bsplines.jl","Residuals.jl","Asymptotics.jl","Norms.jl","Data.jl",
           "Response_moments.jl","Profile_integrals.jl","Profile_response_bounds.jl","Table1_fields.jl",
           "Files.jl","Run_columns.jl","Assemble.jl","Project.toml","Manifest.toml")
    signature=Dict(name=>file_hash(joinpath(@__DIR__,name)) for name in names)
    signature["Profile/Interval_data.jld2"]=file_hash(profile_file())
    return signature
end
