# Data.jl : the data of the computation (Subsection "Data and parameters"), read from Interval_data.jld2, the 128-bit interval
# data prepared once by Prepare_data.jl from the npz files of the numerics.
#
# The interior fields Ψ_in, Ω_in (the plain fields Ψ/z, Ω/z, W = 1) and Ψ̃_in, Ω̃_in (normalized) are tensor-product
# B-splines of degree p on the clamped knot vectors tu, tv of S_L, with interval coefficients enclosing the double-double
# numbers c = hi + lo of the numerics; the integration partitions U, V are the breakpoints of the knots together with the
# point L-d, so that every cell lies in S_{L-d} or in S_band. The knots and the partitions are thin intervals, λ̃ is the
# 128-bit enclosure of the decimal (lambda:app).
function load_data(folder)
    setrounding(BigFloat, RoundNearest)                                # JLD2 restores the BigFloat bounds from decimal strings
    data = load(joinpath(folder, "Interval_data.jld2"), "data")
    precision(BigFloat) == data.prec || error("the data are $(data.prec)-bit intervals and the working precision is $(precision(BigFloat)) bits")
    return validate_data(data)
end

function validate_data(data)
    p = data.p; L = data.L; d = data.d; tu = data.tu; tv = data.tv
    p isa Integer && p ≥ 2 && 0 < d < L || error("invalid spline degree or cutoff parameters")
    length(tu) ≥ 2p+2 && length(tv) ≥ 2p+2 || error("not enough knots")
    for x in (tu, tv, data.U, data.V, data.Ψ, data.Ω, data.Ψ̃, data.Ω̃)
        all(isguaranteed, x) && all(y -> isfinite(inf(y)) && isfinite(sup(y)) && precision(inf(y)) == precision(sup(y)) == data.prec, x) || error("the data are not guaranteed finite intervals at $(data.prec) bits")
    end
    all(isthin, tu) && all(isthin, tv) && issorted(inf.(tu)) && issorted(inf.(tv)) || error("the knots are not thin and sorted")
    all(sup.(tu[1:p+1]) .== 0) && all(inf.(tu[end-p:end]) .== L) && all(sup.(tv[1:p+1]) .== 0) && all(inf.(tv[end-p:end]) .== L) || error("the knots are not clamped on [0,L] with L = $L")
    for t in (tu, tv), x in unique(inf.(t))
        (x == 0 || x == L || count(y -> inf(y) == x, t) ≤ p-1) || error("an interior knot does not give a continuously differentiable spline")
    end
    size(data.Ψ) == size(data.Ω) == size(data.Ψ̃) == size(data.Ω̃) == (length(tu) - p - 1, length(tv) - p - 1) || error("the coefficient arrays do not match the knots")
    for (X, t) in ((data.U, tu), (data.V, tv))
        length(X) ≥ 2 && all(isthin, X) && all(diff(inf.(X)) .> 0) && sup(X[1]) == 0 && inf(X[end]) == L && any(x -> inf(x) == L - d, X) || error("the partition is not a thin strictly increasing partition of [0,L] containing L-d")
        all(x -> inf(x) in inf.(X), unique(t)) || error("a partition cell crosses a spline knot")
    end
    isguaranteed(data.λ̃) && isfinite(inf(data.λ̃)) && isfinite(sup(data.λ̃)) && precision(inf(data.λ̃)) == precision(sup(data.λ̃)) == data.prec || error("λ̃ is not a guaranteed finite $(data.prec)-bit interval")
    in_interval(big(235059685902438)//big(10)^15, data.λ̃) || error("λ̃ does not enclose the exact decimal (lambda:app)")
    return (p = p, tu = tu, tv = tv, U = data.U, V = data.V, L = L, d = d, λ̃ = data.λ̃, r★ = sqrt(U0(interval(BigFloat, L - d))),
            Ψ = data.Ψ, Ω = data.Ω, Ψ̃ = data.Ψ̃, Ω̃ = data.Ω̃)
end

