# Taylor.jl : local polynomial models on the normalized cell [-1,1]^2 (Section "Local Taylor expansion").
#
# A model of a function g on the cell is M = Model(Π, ϱ, N), where Π is a polynomial with interval coefficients
# (a Sequence of RadiiPolynomial in Taylor(n1) ⊗ Taylor(n2), coordinates (x,y) of (cap:local-coordinates)),
# ϱ ≥ 0 bounds |g - Π| on the cell, and N is the bidegree at which the products are truncated.
# In the interior N exceeds every degree that occurs, so that the arithmetic is exact (ϱ = 0 throughout).

struct Model
    Π::Sequence
    ϱ::Interval
    N::Int
    function Model(Π::Sequence,ϱ::Interval,N::Int)
        P=eltype(coefficients(Π))==Interval{Float64} ? Π : Sequence(space(Π),interval.(Float64,coefficients(Π)))
        e=interval(Float64,ϱ)
        isguaranteed(e) && isfinite(inf(e)) && isfinite(sup(e)) && inf(e)≥0 || error("invalid Float64 model error")
        new(P,e,N)
    end
end
Model(Π::Sequence, N::Int) = Model(Π, I"0", N)     # exact polynomial

import Base: +, -, *, /, ^, log

+(M1::Model, M2::Model) = Model(M1.Π + M2.Π, M1.ϱ + M2.ϱ, min(M1.N, M2.N))
-(M1::Model, M2::Model) = Model(M1.Π - M2.Π, M1.ϱ + M2.ϱ, min(M1.N, M2.N))
-(M::Model) = Model(-M.Π, M.ϱ, M.N)
*(c::Interval, M::Model) = (d=interval(Float64,c); Model(d*M.Π, abs(d)*M.ϱ, M.N))
*(M::Model, c::Interval) = c*M
*(c::Integer, M::Model) = interval(c)*M
# Form the reciprocal in the coefficient precision, also when c is an I"..." Float64 interval.
/(M::Model, c::Interval) = (one(M.Π[(0,0)])/interval(Float64,c))*M
+(c::Interval, M::Model) = Model(M.Π + interval(Float64,c), M.ϱ, M.N)
+(M::Model, c::Interval) = c + M
-(c::Interval, M::Model) = c + (-M)
-(M::Model, c::Interval) = M + (-c)

# product of two models, (cap:tm-product), truncated to the bidegree (N,N) with N = min(N1,N2)
function *(M1::Model, M2::Model)
    N = min(M1.N, M2.N)
    Π, τ = πN(M1.Π*M2.Π, N)
    return Model(Π, M1.ϱ*norm(M2.Π,1) + M2.ϱ*norm(M1.Π,1) + M1.ϱ*M2.ϱ + τ, N)
end

function ^(M::Model, n::Integer)
    n ≥ 1 || error("integer powers of models: the exponent must be positive")
    P = M
    for k = 2:n
        P = P*M
    end
    return P
end

# π_N Π and |(I-π_N)Π|_1, (cap:coefficient-norm), with π_N the projection on the bidegree ≤ (N,N): the sum of the
# absolute values of the discarded coefficients only (the difference Π - π_N Π would add the widths of the kept ones)
function πN(Π, N)
    n1, n2 = order(space(Π))
    if n1 ≤ N && n2 ≤ N
        return Π, I"0"
    end
    ΠN = project(Π, Taylor(min(n1,N)) ⊗ Taylor(min(n2,N)))
    τ = I"0"
    for (i,j) in indices(space(Π))
        if i > N || j > N
            τ = τ + abs(Π[(i,j)])
        end
    end
    return ΠN, τ
end

# sup |g| ≤ |Π|_1 + ϱ on the cell, Lemma l:tmops
sup_bound(M::Model) = norm(M.Π,1) + M.ϱ

# model of g^2 without truncation (used before integrating)
square(M::Model) = Model(M.Π^2, I"2"*M.ϱ*norm(M.Π,1) + M.ϱ*M.ϱ, M.N)

# model of g^α for a model M of g > 0, Lemma l:binom, with c the constant coefficient of Π, η = (g-c)/c,
# S_N = c^α Σ_{k≤N} binom(α,k) η^k computed by Horner and the tail (cap:tm-power)
function ^(M::Model, α::Interval)
    N = M.N
    c = interval(mid(M.Π[(0,0)]))
    inf(c) > 0 || error("binomial series: the constant coefficient is not positive")
    α = one(c)*α                                                  # promote exact I"..." exponents before forming a noninteger ratio
    η = (I"1"/c)*(M - c)
    q = norm(η.Π,1) + η.ϱ
    ρ = q*max(I"1", abs(α - interval(N))/interval(N+1))
    sup(ρ) < 1 || error("binomial series: ρ_N ≥ 1, the cell must be subdivided")
    binom = [one(c)]                                               # binom(α,k) = binom(α,k-1) (α-k+1)/k, in the precision of c
    for k = 1:N+1
        push!(binom, binom[k]*(α - interval(k-1))/interval(k))
    end
    S = Model(binom[N+1]*one(M.Π), N)
    for k = N-1:-1:0
        S = binom[k+1] + η*S
    end
    tail = c^α*abs(binom[N+2])*q^interval(N+1)/(I"1" - ρ)
    return Model(c^α*S.Π, c^α*S.ϱ + tail, N)
end

# model of log g for a model M of g > 0, Lemma l:log, with S = log c + Σ_{k=1}^N (-1)^{k+1} η^k/k and the tail (C:logeps)
function log(M::Model)
    N = M.N
    c = interval(mid(M.Π[(0,0)]))
    inf(c) > 0 || error("logarithm: the constant coefficient is not positive")
    η = (I"1"/c)*(M - c)
    q = norm(η.Π,1) + η.ϱ
    sup(q) < 1 || error("logarithm: q ≥ 1, the cell must be subdivided")
    S = Model(one(c)*interval((-1)^(N+1))/interval(N)*one(M.Π), N)   # the coefficients (-1)^{k+1}/k in the precision of c
    for k = N-1:-1:1
        S = one(c)*interval((-1)^(k+1))/interval(k) + η*S
    end
    S = log(c) + η*S
    return Model(S.Π, S.ϱ + q^interval(N+1)/(interval(N+1)*(I"1" - q)), N)
end

# rigorous conversion of a model to Float64 intervals
float64(M::Model) = Model(Sequence(space(M.Π), interval.(Float64, coefficients(M.Π))), interval(Float64, M.ϱ), M.N)

# ∫_{[-1,1]^2} Π dx dy for a polynomial Π, (cap:tm-integral)
function integral(Π)
    I = I"0"
    for (i,j) in indices(space(Π))
        if iseven(i) && iseven(j)
            I = I + Π[(i,j)]*I"4"/interval(i+1)/interval(j+1)
        end
    end
    return I
end

# the moments m_i = ∫_{-1}^1 x^i Π(x) dx, i = 0,...,d, of a polynomial Π(x) stored in Taylor(n) ⊗ Taylor(0)
function moments(Π, d)
    n = order(space(Π))[1]
    m = zeros(eltype(coefficients(Π)), d+1)
    for i = 0:d, k = 0:n
        if iseven(i+k)
            m[i+1] = m[i+1] + Π[(k,0)]*I"2"/interval(i+k+1)
        end
    end
    return m
end

# a density ρ(x) ≥ 0 on the cell, given by a model, prepared for the integrals: its moments, |Π_ρ|_1 and ϱ_ρ.
# The polynomial of a density in y is stored in Taylor(0) ⊗ Taylor(n): it is transposed to Taylor(n) ⊗ Taylor(0) first.
function density(M::Model, d)
    n1, n2 = order(space(M.Π))
    n1 == 0 || n2 == 0 || error("a density is a function of one variable")
    Π = M.Π
    if n2 > 0
        Π = zeros(eltype(coefficients(M.Π)), Taylor(n2) ⊗ Taylor(0))
        for k = 0:n2
            Π[(k,0)] = M.Π[(0,k)]
        end
    end
    return (m = moments(Π, d), Π1 = norm(Π,1), ϱ = M.ϱ)
end

# enclosure of ∫_{[-1,1]^2} g(x,y)^2 ρ(x) μ(y) dx dy for a model M of g and the densities ρ, μ ≥ 0 prepared by
# density(), through the model Q = (Π_Q, ϱ_Q) = square(M) of g^2, whose polynomial part Π_Q = Π_M^2 is ≥ 0 on the cell:
# the polynomial part is integrated exactly by the moments, and the function parts are bounded by ϱ_{ρμ} ∫ Π_Q and
# ϱ_Q ∫ ρμ ≤ ϱ_Q (m_0 n_0 + 4 ϱ_{ρμ}) (Lemma cap:integral-lemma in the form (e:int) of the companion note)
function cell_integral(M::Model, ρ, μ)
    Q = square(M)
    I = I"0"
    for (i,j) in indices(space(Q.Π))
        I = I + Q.Π[(i,j)]*ρ.m[i+1]*μ.m[j+1]
    end
    ϱ = ρ.ϱ*μ.Π1 + μ.ϱ*ρ.Π1 + ρ.ϱ*μ.ϱ
    I = I + I"[-1,1]"*(ϱ*integral(Q.Π) + Q.ϱ*(ρ.m[1]*μ.m[1] + I"4"*ϱ))
    return intersect_interval(I, I"[0,inf]")
end

# the polynomial Π(x,y(σ)) with y = (σ^2 + 2σ - 1)/2, i.e. the chart (cap:axis-chart) v = s^2 with s = (√d/2)(1+σ),
# on a cell touching the axis v = 0 (the model's remainder is unchanged by the substitution)
function substitute_axis(M::Model)
    n1, n2 = order(space(M.Π))
    T = eltype(coefficients(M.Π))
    y = zeros(T, Taylor(0) ⊗ Taylor(2))
    y[(0,0)] = I"-0.5"; y[(0,1)] = I"1"; y[(0,2)] = I"0.5"
    Π = [zeros(T, Taylor(n1) ⊗ Taylor(0)) for j = 0:n2]            # Π(x,y) = Σ_j Π_j(x) y^j
    for j = 0:n2, i = 0:n1
        Π[j+1][(i,0)] = M.Π[(i,j)]
    end
    S = Π[n2+1]
    for j = n2-1:-1:0                                             # Horner in y
        S = S*y + Π[j+1]
    end
    return Model(S, M.ϱ, M.N)
end
