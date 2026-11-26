# Input files keep their original exact knots and certified spectral data.
# Every coefficient used by a local Model or a matrix is rounded outward to Float64.
proof_float64(x::Interval) = interval(Float64,x)
proof_float64(x::Sequence) = Sequence(space(x),proof_float64.(coefficients(x)))
proof_float64(x::NamedTuple) = NamedTuple{keys(x)}(map(proof_float64,values(x)))
proof_float64(x::AbstractArray) = map(proof_float64,x)
proof_float64(x::Tuple) = map(proof_float64,x)
proof_float64(x::AbstractDict) = Dict(k=>proof_float64(v) for (k,v) in x)
proof_float64(x) = x
