include(joinpath(@__DIR__,"..","Paths.jl"))
# Use the same interval arithmetic and local polynomials as the Table 2 computation.
using IntervalArithmetic, RadiiPolynomial, LinearAlgebra, JLD2, SHA
setprecision(BigFloat,128)
setrounding(BigFloat,RoundNearest)
IntervalArithmetic.configure(rounding=:correct,matmul=:slow)
BLAS.set_num_threads(1)

include("Float64.jl")

const CORE = joinpath(@__DIR__,"reused")
const OUTPUT = joinpath(COERCIVITY_RESULTS,"proposition")
for file in ("Taylor.jl","Bsplines.jl","Residuals.jl","Asymptotics.jl","Norms.jl","Data.jl",
             "Functionals.jl","Response_data.jl","Response_fields.jl","Response_moments.jl")
    include(joinpath(CORE,file))
end

isdefined(@__MODULE__,:load_coercivity_proposition) || include(joinpath(@__DIR__,"..","Compact_data.jl"))

include("Profile.jl")
include("Boxes.jl")
include("Positivity.jl")
include("Farfield.jl")
include("Matrices.jl")
include("Data.jl")
include("Spectral.jl")
include("Gram.jl")
include("Grids.jl")
include("Bound_columns.jl")
include("Matrix_columns.jl")
include("Gram_columns.jl")
