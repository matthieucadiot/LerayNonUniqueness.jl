# Bsplines.jl : local polynomials of the tensor-product B-spline representation (Section "Spline representations").
#
# Thin 128-bit knots and cell endpoints fix the topology, including tiny junction cells.
# Local spline conversion uses 256 bits to preserve its cancellation; the finished
# polynomial coefficients are rounded outward once to Float64 before model arithmetic.
# The knots, the cell endpoints and the coefficients are the intervals of Interval_data.jld2 (Data.jl).

# The degree-p basis functions B_{k-p,p},...,B_{k,p} on the knot interval [t_k,t_{k+1}] ⊇ [a,b], as polynomials of
# the normalized variable x ∈ [-1,1], u = (a+b)/2 + (b-a)/2 x, by the recursion (eq:coxdeboor). The knots t and the
# endpoints a, b are thin intervals, the result is the (p+1)×(p+1) matrix B with B[n,i+1] the coefficient of x^i of B_{k-p+n-1,p}.
function local_basis(t, p, k, a, b)
    return setprecision(2*precision(BigFloat)) do
        u = Sequence(Taylor(1), [(a + b)/I"2", (b - a)/I"2"])
        Bq = [Sequence(Taylor(0), [big(I"1")])]                 # B_{k,0} = 1 on [t_k,t_{k+1})
        for q = 1:p
            Bnew = Vector{Any}(undef, q+1)
            for m = 1:q+1
                n = k - q + m - 1                                       # B_{n,q} from B_{n,q-1} = Bq[m-1] and B_{n+1,q-1} = Bq[m]
                P = zeros(Interval{BigFloat}, Taylor(q))
                if m ≥ 2 && inf(t[n+q]) > sup(t[n])                     # the repeated knots give no term
                    P = P + (u - t[n])*Bq[m-1]/(t[n+q] - t[n])
                end
                if m ≤ q && inf(t[n+q+1]) > sup(t[n+1])
                    P = P + (t[n+q+1] - u)*Bq[m]/(t[n+q+1] - t[n+1])
                end
                Bnew[m] = P
            end
            Bq = Bnew
        end
        [Bq[n][i] for n = 1:p+1, i = 0:p]
    end
end

# index k of the knot interval [t_k,t_{k+1}) containing the point a (the last knot ≤ a), for thin knots and a thin a
knot_interval(t, a) = searchsortedlast(inf.(t), inf(a))

# the polynomial Σ c_{n1,n2} B_{n1}(u) B_{n2}(v) of a spline on the cell, (eq:tensorloc), in the normalized variables,
# with Bu, Bv the local bases of local_basis, ku, kv the knot interval indices, and c the interval coefficients of the spline
function piece(c, Bu, ku, Bv, kv, p, denominator=1)
    prec = precision(BigFloat)
    M = setprecision(2*prec) do
        A = transpose(Bu)*c[ku-p:ku, kv-p:kv]*Bv
        # A common exact denominator keeps the saved coefficient numerators
        # intact throughout the de Boor basis conversion on microscopic cells.
        denominator == 1 ? A : A ./ interval(BigFloat,denominator)
    end
    all(isguaranteed, M) || error("the spline piece is not guaranteed")
    Π = zeros(Interval{Float64}, Taylor(p) ⊗ Taylor(p))
    for i = 0:p, j = 0:p
        Π[(i,j)] = interval(Float64,M[i+1,j+1])                                  # explicit outward rounding of the fixed local polynomial
    end
    return Π
end

# the field and its derivatives ∂_u, ∂_v, ∂_u^2, ∂_v^2 on the cell as exact models, with ∂_u = ∂_x/h_u, ∂_v = ∂_y/h_v
function field(Π, hu, hv, N)
    Πu = differentiate(Π, (1,0))/hu
    Πv = differentiate(Π, (0,1))/hv
    Πuu = differentiate(Πu, (1,0))/hu
    Πvv = differentiate(Πv, (0,1))/hv
    return (f = Model(Π, N), u = Model(Πu, N), v = Model(Πv, N), uu = Model(Πuu, N), vv = Model(Πvv, N))
end
