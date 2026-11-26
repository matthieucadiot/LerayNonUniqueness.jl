# Asymptotics.jl : the asymptotic solutions and the exterior (Subsections "Arithmetic at infinity" and
# "Computation on an exterior domain").
#
# A far-field expression is a polynomial F(x,t,ℓ) = Σ F_{khj} x^k t^h ℓ^j of (cap:asymptotic-polynomial), stored as
# a Sequence in Taylor(K) ⊗ Taylor(H) ⊗ Taylor(J), with x = r^2/R^2, t = 1/R, ℓ = log R, so that
# r^2 = x R^2, z^2 = (1-x) R^2, and the homogeneous angular factor R^2/R^2 is represented by 1.

# the variables x, t, ℓ as polynomials (BigFloat intervals at the current precision)
function variables(endpoint=Float64)
    X = zeros(Interval{endpoint}, Taylor(1) ⊗ Taylor(0) ⊗ Taylor(0)); X[(1,0,0)] = I"1"
    T = zeros(Interval{endpoint}, Taylor(0) ⊗ Taylor(1) ⊗ Taylor(0)); T[(0,1,0)] = I"1"
    ℓ = zeros(Interval{endpoint}, Taylor(0) ⊗ Taylor(0) ⊗ Taylor(1)); ℓ[(0,0,1)] = I"1"
    return X, T, ℓ
end

# 𝒜_r = 2(1-x)∂_x - t∂_t + ∂_ℓ and 𝒜_z = -2x∂_x - t∂_t + ∂_ℓ, (cap:asymptotic-derivatives): (∂_r F)/r = t^2 𝒜_r F, (∂_z F)/z = t^2 𝒜_z F
function 𝒜r(F)
    X, T, ℓ = variables(typeof(inf(coefficients(F)[1])))
    return I"2"*(I"1" - X)*differentiate(F, (1,0,0)) - T*differentiate(F, (0,1,0)) + differentiate(F, (0,0,1))
end
function 𝒜z(F)
    X, T, ℓ = variables(typeof(inf(coefficients(F)[1])))
    return -I"2"*X*differentiate(F, (1,0,0)) - T*differentiate(F, (0,1,0)) + differentiate(F, (0,0,1))
end

# 𝓛_ε F = t^2 ((4-2x) 𝒜_r F + (2x-1+2ε) 𝒜_z F + x 𝒜_r^2 F + (1-x) 𝒜_z^2 F), (cap:asymptotic-laplacian):
# L_0 F = 𝓛_0 F for an even field, L_0(zF)/z = 𝓛_1 F for an odd field zF
function 𝓛(F, ε)
    ε in (0, 1) || error("far-field Laplacian: ε must be 0 or 1")
    X, T, ℓ = variables(typeof(inf(coefficients(F)[1])))
    ArF = 𝒜r(F)
    AzF = 𝒜z(F)
    return T*T*((I"4" - I"2"*X)*ArF + (I"2"*X - I"1" + I"2"*interval(ε))*AzF + X*𝒜r(ArF) + (I"1" - X)*𝒜z(AzF))
end

is_thin(F) = all(inf(c) == sup(c) for c in coefficients(F))
is_zero(F) = all(inf(c) == 0 && sup(c) == 0 for c in coefficients(F))
contains_zero(F) = all(inf(c) ≤ 0 && sup(c) ≥ 0 for c in coefficients(F))

# the far-field polynomial with its bounds rounded outward to the working precision: the numbers built in 1024-bit
# arithmetic keep their precision through every product (IntervalArithmetic rounds at the precision of the operands),
# which would make the arithmetic on the band cells needlessly slow
function working_precision(F)
    all(isguaranteed, coefficients(F)) || error("far field: coefficients are not guaranteed intervals")
    prec = precision(BigFloat)
    return Sequence(space(F),interval.(Float64,coefficients(F)))
end

# The asymptotic solutions (e:trnl), (e:trnl2) at the constants (constants), (constants2), divided by z for the odd
# fields: P = Ψ_{≥-9}/z, O = Ω_{≥-11}/z, P̃ = Ψ̃_{≥-7}, Õ = Ω̃_{≥-9}, with the formulas (gf2)-(gf29), and the exterior
# residual E_Ω/z of (sec3eq7.6), i.e. the second equation of (gf1) at Ψ = zP, Ω = zO with r∂_r = x 𝒜_r, z∂_z = (1-x) 𝒜_z:
#   E_Ω/z = 𝓛_1 O + (x 𝒜_r O + (1-x) 𝒜_z O + 4 O)/2 + (x 𝒜_r P + 2P)(O + (1-x) 𝒜_z O) - x (P + (1-x) 𝒜_z P) 𝒜_r O.
# The odd fields are built times K = 213050462625, a common multiple of the denominators of their rational
# coefficients as written below (σ_2 = 1440081/5 included), in 1024-bit arithmetic, where every operation below is an exact integer
# operation: the identity 𝓛_1 P = O and the vanishing of the orders t^h, h ≤ 12, of E_Ω/z are then obtained exactly.
# Final coefficient bounds are rounded outward to the working precision prec.
function far_field(prec)
    return setprecision(1024) do
    X, T, ℓ = variables(BigFloat)
    r2 = X; z2 = I"1" - X; R2 = one(X)                     # r^2/R^2, z^2/R^2, R^2/R^2
    σ = big(I"130"); σ₂ = big(I"1440081"); σ₃ = big(I"90012367"); σ₄ = big(I"32103000000")   # σ_2 = σ₂/5
    σ² = σ*σ; σ³ = σ²*σ; σ⁴ = σ³*σ; σ⁵ = σ⁴*σ
    h4 = r2^2 - I"12"*r2*z2 + I"8"*z2^2                             # h_4/R^4, (gf25)
    q4 = I"5"*r2^2 - I"20"*r2*z2 + I"8"*z2^2                 # (5r^4-20r^2z^2+8z^4)/R^4

    # K Ψ_{≥-9}/z: each line is K n/d times a polynomial with integer coefficients, (gf2), (gf6), (gf7), (gf10), (gf13), (gf17), (gf20), (gf23)
    KP = I"213050462625"*σ*T^2 +
         I"71016820875"*T^4*(σ²*(I"3"*r2 - I"7"*z2) - I"18"*σ*R2) +
         I"42610092525"*σ₂*T^5 +
         I"71016820875"*T^6*(σ³*(-I"8"*r2^2 - I"47"*r2*z2 + I"29"*z2^2) - I"108"*σ²*R2*(r2 - z2) + I"36"*σ*R2^2) +
         I"405810405"*T^7*(I"28"*σ*σ₂*ℓ*(I"3"*r2 - I"4"*z2) + I"525"*σ₃*(I"3"*r2 - I"4"*z2) + I"36"*σ*σ₂*R2) +
         I"14203364175"*T^8*(σ⁴*(-I"451"*r2^3 - I"149"*r2^2*z2 + I"2685"*r2*z2^2 - I"537"*z2^3)
                            + I"6"*σ³*R2*(-I"471"*r2^2 + I"1462"*r2*z2 - I"387"*z2^2)
                            + I"36"*σ²*R2^2*(I"27"*r2 - I"41"*z2) + I"360"*σ*R2^3) +
         T^9*(R2^2*(I"3298014720"*σ²*σ₂*ℓ - I"36011097980"*σ²*σ₂ + I"37566448920"*σ*σ₂ + I"61837776000"*σ*σ₃)
              + R2*(I"3"*r2 - I"4"*z2)*(I"2529727200"*σ²*σ₂*ℓ + I"32716436220"*σ²*σ₂ - I"4553508960"*σ*σ₂ + I"47432385000"*σ*σ₃)
              + q4*(I"1748106360"*σ²*σ₂*ℓ^2 + I"29282964480"*σ²*σ₂*ℓ + I"65553988500"*σ*σ₃*ℓ + I"213050462625"*σ₄)) +
         I"27054027"*T^10*(σ⁵*(-I"1413624"*r2^4 + I"5108254"*r2^3*z2 + I"9506631"*r2^2*z2^2 - I"11069496"*r2*z2^3 + I"785751"*z2^4)
                               + I"36"*σ⁴*R2*(-I"202799"*r2^3 + I"1280103"*r2^2*z2 - I"1106397"*r2*z2^2 + I"109201"*z2^3)
                               + I"36"*σ³*R2^2*(I"43578"*r2^2 - I"224719"*r2*z2 + I"117203"*z2^2)
                               + I"151200"*σ²*R2^3*(I"9"*r2 - I"16"*z2) + I"1134000"*σ*R2^4)

    # K Ω_{≥-11}/z, (gf2), (gf4), (gf9), (gf12), (gf15), (gf19), (gf22)
    KO = -I"6"*I"213050462625"*σ*T^4 +
         I"213050462625"*T^6*(I"12"*σ²*(-I"2"*r2 + I"3"*z2) + I"24"*σ*R2) +
         I"213050462625"*T^8*(I"2"*σ³*(-I"39"*r2^2 + I"190"*r2*z2 - I"111"*z2^2) + I"72"*σ²*R2*(I"3"*r2 - I"5"*z2) + I"72"*σ*R2^2) +
         I"102264222060"*σ*σ₂*T^9*(-r2 + I"6"*z2) +
         I"213050462625"*T^10*(I"8"*σ⁴*(-I"30"*r2^3 + I"282"*r2^2*z2 - I"531"*r2*z2^2 + I"179"*z2^3)
                          + I"16"*σ³*R2*(I"78"*r2^2 - I"403"*r2*z2 + I"215"*z2^2)
                          + I"96"*σ²*R2^2*(I"6"*r2 - I"11"*z2) + I"576"*σ*R2^3) +
         T^11*(I"19478899440"*σ²*σ₂*ℓ*(I"3"*r2^2 + I"56"*r2*z2 - I"24"*z2^2)
               + I"365229364500"*σ*σ₃*(I"3"*r2^2 + I"56"*r2*z2 - I"24"*z2^2)
               - I"4174049880"*σ²*σ₂*(I"279"*r2^2 - I"1001"*r2*z2 + I"1702"*z2^2)
               + I"350620189920"*σ*σ₂*R2*(I"3"*r2 + I"5"*z2)) +
         I"42610092525"*T^12*(I"2"*σ⁵*(-I"1938"*r2^4 + I"24458"*r2^3*z2 - I"99933"*r2^2*z2^2 + I"104448"*r2*z2^3 - I"23623"*z2^4)
                            + I"8"*σ⁴*R2*(I"3403"*r2^3 - I"39831"*r2^2*z2 + I"62949"*r2*z2^2 - I"17177"*z2^3)
                            + I"24"*σ³*R2^2*(I"293"*r2^2 - I"2764"*r2*z2 + I"2083"*z2^2)
                            + I"28800"*σ²*R2^3*(r2 - I"2"*z2) + I"36000"*σ*R2^4)

    is_thin(KP) && is_thin(KO) || error("far field: the arithmetic is not exact")
    is_zero(𝓛(KP, 1) - KO) || error("far field: 𝓛_1 P ≠ O")

    # 2K² E_Ω/z, exact; the orders t^h, h ≤ 12, vanish identically
    ArP = 𝒜r(KP); AzP = 𝒜z(KP); ArO = 𝒜r(KO); AzO = 𝒜z(KO)
    E = I"213050462625"*(I"2"*𝓛(KO, 1) + X*ArO + z2*AzO + I"4"*KO) + I"2"*(X*ArP + I"2"*KP)*(KO + z2*AzO) - I"2"*X*(KP + z2*AzP)*ArO
    is_thin(E) || error("far field: the arithmetic of the exterior residual is not exact")
    for (k,h,j) in indices(space(E))
        if h ≤ 12 && !(inf(E[(k,h,j)]) == 0 && sup(E[(k,h,j)]) == 0)
            error("far field: the exterior residual has a term of order t^$h")
        end
    end

    # the even fields (gf26), (gf29), (gf28) at the constants (constants2), with variables and constants built at the
    # working precision, and the odd fields divided by K and rounded outward to the working precision
    X, T, ℓ = variables(BigFloat)
    r2 = X; z2 = I"1" - X
    σ = big(I"130"); h4 = r2^2 - I"12"*r2*z2 + I"8"*z2^2
    τ₂ = big(I"22794186")/big(I"100000"); τ₃ = -big(I"15237197")/big(I"100"); λ̃ₐ = big(I"235059686")/big(I"1000000000")
    P̃ = τ₂*T^5*(I"4"*z2 - r2) - σ*τ₂/(λ̃ₐ + I"3")*T^7*(big(I"4")/big(I"49")*(I"17"*r2^2 - I"15"*r2*z2 - I"32"*z2^2) + big(I"90")/big(I"77")*ℓ*h4) + τ₃*T^7*h4
    Õ = -I"18"*σ*τ₂/(λ̃ₐ + I"3")*T^9*(r2^2 + I"8"*r2*z2 - I"8"*z2^2)
    contains_zero(𝓛(P̃, 0) - Õ) || error("far field: 𝓛_0 P̃ ≠ Õ")
    return (P = working_precision(KP/big(I"213050462625")), O = working_precision(KO/big(I"213050462625")), P̃ = working_precision(P̃), Õ = working_precision(Õ),
            EΩ = working_precision(E/(big(I"2")*I"213050462625"*I"213050462625")))
    end
end

# E_Ω̃ of (sec3eq7.6) in the exterior: the second equation of (gf24) at the eigenvalue λ̃, with Ψ = zP, Ω = zO, Ψ̃ = P̃, Ω̃ = Õ,
#   E_Ω̃ = 𝓛_0 Õ + (x 𝒜_r Õ + (1-x) 𝒜_z Õ + (3-2λ̃) Õ)/2 - x (P + (1-x) 𝒜_z P) 𝒜_r Õ + (1-x)(x 𝒜_r P + 2P) 𝒜_z Õ
#         - x(1-x) 𝒜_z P̃ 𝒜_r O + (x 𝒜_r P̃ + 2P̃)(O + (1-x) 𝒜_z O)
function exterior_residual_even(far, λ̃)
    X, T, ℓ = variables()
    z2 = I"1" - X
    P = far.P; O = far.O; P̃ = far.P̃; Õ = far.Õ
    return 𝓛(Õ, 0) + I"0.5"*(X*𝒜r(Õ) + z2*𝒜z(Õ) + (I"3" - I"2"*λ̃)*Õ) - X*(P + z2*𝒜z(P))*𝒜r(Õ) + z2*(X*𝒜r(P) + I"2"*P)*𝒜z(Õ) - X*z2*𝒜z(P̃)*𝒜r(O) + (X*𝒜r(P̃) + I"2"*P̃)*(O + z2*𝒜z(O))
end

# the far-field polynomials of F, (∂_r F)/r = t^2 𝒜_r F, (∂_r/r)^2 F = t^4 (𝒜_r^2 F - 2 𝒜_r F), and the same in z
# (Proposition p:ops), used by far_field_on_cell in the gluing region
function far_field_derivatives(F)
    F=working_precision(F)
    X, T, ℓ = variables()
    Ar = 𝒜r(F)
    Az = 𝒜z(F)
    return (f = F, r1 = T*T*Ar, r2 = T^4*(𝒜r(Ar) - I"2"*Ar), z1 = T*T*Az, z2 = T^4*(𝒜z(Az) - I"2"*Az))
end

# ---------------------------------------------------------------------------------------------------------------
# The exterior integrals and suprema, at the active BigFloat precision

# rigorous conversion of a far-field polynomial to Float64 intervals
float64(F::Sequence) = Sequence(space(F), interval.(Float64, coefficients(F)))

# 𝒜_{m,n} = ∫_0^{π/2} sin^m θ cos^n θ dθ by the recursions following (X:quadform)
function 𝒜(m::Integer, n::Integer)
    m ≥ 0 && n ≥ 0 || error("angular moment: m and n must be nonnegative")
    if m ≥ 2
        return interval(Float64, m-1)/interval(Float64, m+n)*𝒜(m-2, n)
    elseif n ≥ 2
        return interval(Float64, n-1)/interval(Float64, m+n)*𝒜(m, n-2)
    elseif m == 1 && n == 1
        return I"0.5"
    elseif m == 0 && n == 0
        return interval(Float64, π)/I"2"
    else
        return I"1"
    end
end

# ℬ_{q,j}(L) = ∫_L^∞ R^q (log R)^j dR for q < -1, (cap:exterior-moments)
function ℬ(q::Integer, j::Integer, L)
    L = interval(Float64, L)
    q < -1 && j ≥ 0 && inf(L) > 0 || error("radial moment: q < -1, j ≥ 0 and L > 0 are required")
    S = I"0"
    for i = 0:j
        factor = factorial(big(j))÷factorial(big(j-i))                   # exact integer, avoiding machine-integer overflow
        S = S + interval(Float64, factor)*log(L)^interval(Float64, j-i)/interval(Float64, -q-1)^interval(Float64, i)
    end
    return L^interval(Float64, q+1)/interval(Float64, -q-1)*S
end

# 4π² ∫∫_{r,z>0, R≥L} r^{3+a} z^b F^2 dr dz for a far-field polynomial F, Proposition p:gram; the coefficients of F^2
# keep their signs, and every retained term is checked to be integrable
function exterior_integral(F, a::Integer, b::Integer, L)
    F=working_precision(F); L=interval(Float64,L)
    a ≥ 0 && b ≥ 0 && inf(L) > 0 || error("exterior integral: a, b ≥ 0 and L > 0 are required")
    F2 = F^2
    I = I"0"
    for (k,h,j) in indices(space(F2))
        c = F2[(k,h,j)]
        if !(inf(c) == 0 && sup(c) == 0)
            h > 5 + a + b || error("exterior integral: the term of order t^$h is not integrable")
            I = I + c*𝒜(2*k+3+a, b)*ℬ(4+a+b-h, j, L)
        end
    end
    return I"4"*interval(Float64, π)*interval(Float64, π)*I
end

# sup_{R≥L} |r z^b F| ≤ Σ |F_{khj}| max_θ(sin^{2k+1}θ cos^b θ) L^{1+b-h} (log L)^j, (cap:exterior-supremum), L > 1
function exterior_sup(F, b::Integer, L)
    F=working_precision(F)
    L = interval(Float64, L)
    b ≥ 0 || error("exterior supremum: b must be nonnegative")
    inf(L) > 1 || error("exterior supremum: L must be larger than 1")
    S = I"0"
    for (k,h,j) in indices(space(F))
        c = F[(k,h,j)]
        if !(inf(c) == 0 && sup(c) == 0)
            (1+b-h ≤ 0 && sup(interval(Float64, 1+b-h)*log(L) + interval(Float64, j)) ≤ 0) || error("exterior supremum: the term of order t^$h (log R)^$j is not decreasing")
            m = 2*k+1
            θ = (b == 0) ? I"1" : (interval(Float64, m)/interval(Float64, m+b))^(interval(Float64, m)/I"2")*(interval(Float64, b)/interval(Float64, m+b))^(interval(Float64, b)/I"2")
            S = S + abs(c)*θ*L^interval(Float64, 1+b-h)*log(L)^interval(Float64, j)
        end
    end
    return S
end

# the far field for the computation: the four fields with their derivative forms, the exterior residuals E_Ω/z and
# E_Ω̃ at the eigenvalue λ̃, and the maximal orders (K,H,J) in (x,t,ℓ) of all the expressions
function far_field_data(prec, λ̃)
    far = far_field(prec)
    EΩ̃ = exterior_residual_even(far, λ̃)
    P = far_field_derivatives(far.P); O = far_field_derivatives(far.O); P̃ = far_field_derivatives(far.P̃); Õ = far_field_derivatives(far.Õ)
    K = 0; H = 0; J = 0
    for F in (P, O, P̃, Õ), G in F
        K = max(K, order(space(G))[1]); H = max(H, order(space(G))[2]); J = max(J, order(space(G))[3])
    end
    return (P = P, O = O, P̃ = P̃, Õ = Õ, K = K, H = H, J = J, EΩ = far.EΩ, EΩ̃ = EΩ̃)
end
