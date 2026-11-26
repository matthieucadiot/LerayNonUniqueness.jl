#### list_functions.jl : the functions of main_proof.jl and approximate_eigenpairs.jl (Proposition A.3 of draft5, Appendix A
#### "Computation of the eigenfunctions") : the series pieces of (A.2), their solution operators (Section A.3), the plateau
#### functions, L(λ) of (A.7), the enclosure of polynomials, the number of zeros and the integral 𝒩 (Section A.4).
#### "companion (k)" is the equation (k) of Eigenpairs_RadiiPolynomial.tex, which details the formulas not written in draft5.
#### A series piece is the tuple (c, h, σ, R, q, G) : x = c + hy, R(y) = ρ(c + hy), σ of (A.3), (A.5), q_j = q[j+1], G_s ≥ G[s+1].



############ Section A.3 : the four roots t of Ŝ(t) = S(t)/t⁵ = 126 - 420t + 540t² - 315t³ + 70t⁴ (companion Section 6) ############

function roots_Shat()
    #### Newton from the double precision roots (numerics), then Lemma A.1 in ℂ (norm |·|) on the disc |t - t̄| ≤ R = 10⁻¹⁰
    Ŝ(t) = (((interval(70)*t - interval(315))*t + interval(540))*t - interval(420))*t + interval(126)
    dŜ(t) = ((interval(280)*t - interval(945))*t + interval(1080))*t - interval(420)
    t̄ = [complex(BigFloat(0.9749663273111581), BigFloat(0.37224858264531774)),
         complex(BigFloat(0.9749663273111581), -BigFloat(0.37224858264531774)),
         complex(BigFloat(1.2750336726888385), BigFloat(0.16427693960281098)),
         complex(BigFloat(1.2750336726888385), -BigFloat(0.16427693960281098))]
    for i = 1:4
        for p = 1:6                                                              #### Newton, through the midpoints
            T = complex(interval(real(t̄[i])), interval(imag(t̄[i]))) ; d = dŜ(T)
            s = Ŝ(T)*conj(d)/(real(d)^2 + imag(d)^2)                            #### Ŝ(t̄)/Ŝ'(t̄), 1/z = conj(z)/|z|²
            t̄[i] = t̄[i] - complex(mid(real(s)), mid(imag(s)))
        end
    end
    R = parse(Interval{BigFloat}, "1e-10")
    t = [complex(interval(real(t̄[i])), interval(imag(t̄[i]))) for i = 1:4]
    for i = 1:4
        T = t[i] ; d = dŜ(T)
        A = conj(d)/(real(d)^2 + imag(d)^2) ; A = complex(interval(mid(real(A))), interval(mid(imag(A))))   #### A ≈ 1/Ŝ'(t̄)
        D = complex(real(T) + interval(-sup(R), sup(R)), imag(T) + interval(-sup(R), sup(R)))             #### a box containing the disc
        z = A*Ŝ(T) ; Y = sqrt(real(z)^2 + imag(z)^2)                                                     #### Y ≥ |A Ŝ(t̄)|
        z = interval(1) - A*dŜ(D) ; Z = sqrt(real(z)^2 + imag(z)^2)                                      #### Z ≥ |1 - A Ŝ'(t)| on the disc
        sup(Z) < 1 || error("The condition Z < 1 for the root $i of Ŝ is not verified : the proof is not conclusive.")
        r = sup(Y/(interval(1) - Z))
        r <= inf(R) || error("The condition r ≤ R for the root $i of Ŝ is not verified : the proof is not conclusive.")
        t[i] = complex(real(T) + interval(-r, r), imag(T) + interval(-r, r))
    end
    for i = 1:4
        for j = i+1:4                                                            #### disjoint discs : these are the four roots
            z = t[i] - t[j]
            inf(sqrt(real(z)^2 + imag(z)^2)) > 2*sup(R) || error("The discs of the roots of Ŝ intersect : the proof is not conclusive.")
        end
    end
    return t
end



############ Sections A.1 and A.3 : the series pieces φa, ψa, ψb, φb of (A.2), their polynomials R, q and G of (A.8) ############

function series_pieces(a, x1, x2, x3, x4, b, radial, t, ν, N)
    #### R(y) = ρ(c + hy) = X³τ⁵Ŝ(τ) (X³ dropped for ρ = P), τ = (x - a)/ℓa on the left ramp, τ = (b - x)/ℓb on the right ramp
    ℓa = x2 - a ; ℓb = b - x3
    e0 = Sequence(Taylor(0), [interval(BigFloat, 1)])                           #### the polynomial 1
    y = Sequence(Taylor(1), [interval(BigFloat, 0), interval(BigFloat, 1)])     #### the polynomial y
    Ŝ(τ) = interval(126)*e0 - interval(420)*τ + interval(540)*τ^2 - interval(315)*τ^3 + interval(70)*τ^4
    τa(X) = (X - a*e0)*inv(ℓa)
    τb(X) = (b*e0 - X)*inv(ℓb)
    function piece(c, h, σ, τ, xzeros)
        #### ζ = (x - c)/h for the zeros x of Q (φa, φb : R = |h|⁵y⁵Q) or of R (ψa, ψb) : q_j = -Σ_ζ ζ^{-j-1}, G_j ≥ Σ_ζ (ν/|ζ|)^j/(|ζ| - ν)
        X = c*e0 + h*y
        R = (radial ? X^3 : e0)*τ(X)^5*Ŝ(τ(X))
        q = fill(interval(BigFloat, 0), N+1) ; G = fill(interval(BigFloat, 0), N+1)
        for z in [(x - c)/h for x in xzeros]
            (in_interval(0, real(R(z))) && in_interval(0, imag(R(z)))) ||
                error("A listed zero of R is not a zero of R (Section A.3) : the proof is not conclusive.")
            absz = interval(inf(sqrt(real(z)^2 + imag(z)^2)))                  #### a lower bound of |ζ|
            inf(absz - ν) > 0 || error("The condition ν < |ζ| of (A.8) is not verified : the proof is not conclusive.")
            w = conj(z)/(real(z)^2 + imag(z)^2) ; wj = w                        #### 1/ζ and its powers
            for j = 0:N
                q[j+1] = q[j+1] - real(wj) ; wj = wj*w
                G[j+1] = G[j+1] + pown(ν/absz, j)/(absz - ν)
            end
        end
        return (c, h, σ, R, q, G)
    end
    z0 = radial ? fill(complex(interval(BigFloat, 0), interval(BigFloat, 0)), 3) : Complex{Interval{BigFloat}}[]   #### x = 0 (triple)
    za = fill(complex(a, interval(BigFloat, 0)), 5) ; zb = fill(complex(b, interval(BigFloat, 0)), 5)             #### x = a, b (order five)
    ta = [a + ℓa*t[i] for i = 1:4] ; tb = [b - ℓb*t[i] for i = 1:4]                                               #### x = a + ℓa t, b - ℓb t
    x̄a = (x1 + x2)/interval(2) ; wa = (x2 - x1)/interval(2) ; x̄b = (x3 + x4)/interval(2) ; wb = (x4 - x3)/interval(2)
    φa = piece(a, x1 - a, 4, τa, [z0; ta])              #### y = (x - a)/δa ∈ [0,1], σ = 4
    ψa = piece(x̄a, wa, -1, τa, [z0; za; ta])            #### y = (x - x̄a)/wa ∈ [-1,1], σ = -1
    ψb = piece(x̄b, wb, -1, τb, [z0; zb; tb])            #### y = (x - x̄b)/wb ∈ [-1,1], σ = -1
    φb = piece(b, x4 - b, 4, τb, [z0; tb])              #### h = -δb : y = (b - x)/δb ∈ [0,1], σ = 4
    return φa, ψa, ψb, φb
end



############ Section A.3 : the solution operators U* = 𝒮^u(1,λ), 𝒞^u(1,0,λ), 𝒞^u(0,1,λ) and ∂λU*, for all λ ∈ Λ ############

function solution_operator(piece, Λ, u0, u1, ν, N)
    #### Ū = π^{≤N}U* by the forward substitution of (A.3) or (A.5) in interval arithmetic (for every λ ∈ Λ), and Lemma A.2 with this Ū :
    #### Y = Y_∞(Ū), Z_{1,0} = 0, ‖U* - Ū‖_ν ≤ rU ; the same for V = π^{≤N}∂λU* from (A.12), ‖∂λU* - V‖_ν ≤ rV (companion (9))
    c, h, σ, R, q, G = piece
    μ = Λ*h^2
    Z1_1 = ν*G[1]/interval(N+1+σ) + ν^2*abs(μ)/interval((N+1)*(N+1+σ))
    sup(Z1_1) < 1 || error("The condition Z_{1,1} < 1 (Lemma A.2) is not verified : the proof is not conclusive.")
    u = fill(interval(BigFloat, 0), N+1) ; u[1] = interval(BigFloat, u0) ; u[2] = interval(BigFloat, u1)
    for n = 2:N                                                                  #### n(n+σ)u_n = -Σ_{k=1}^{n-1} k q_{n-1-k} u_k - μ u_{n-2}
        s = μ*u[n-1]
        for k = 1:n-1
            s = s + interval(k)*q[n-k]*u[k+1]
        end
        u[n+1] = -s/interval(n*(n+σ))
    end
    rU = computation_Y_infinity(u, μ, G, σ, ν, N)/(interval(1) - Z1_1)
    v = fill(interval(BigFloat, 0), N+1)
    for n = 2:N                                                                  #### (A.12) : n(n+σ)v_n = -Σ k q_{n-1-k} v_k - μ v_{n-2} - h² u_{n-2}
        s = μ*v[n-1] + h^2*u[n-1]
        for k = 1:n-1
            s = s + interval(k)*q[n-k]*v[k+1]
        end
        v[n+1] = -s/interval(n*(n+σ))
    end
    Yλ = computation_Y_infinity(v, μ, G, σ, ν, N) + h^2*(pown(ν, N+1)*abs(u[N])/interval((N+1)*(N+1+σ)) +
         pown(ν, N+2)*abs(u[N+1])/interval((N+2)*(N+2+σ))) + h^2*ν^2*rU/interval((N+1)*(N+1+σ))
    rV = Yλ/(interval(1) - Z1_1)
    return Sequence(Taylor(N), u), rU, Sequence(Taylor(N), v), rV
end

function computation_Y_infinity(v, μ, G, σ, ν, N)
    #### Y_∞(v) ≥ ‖π^{>N}f_λ(v)‖_ν of (A.10) for v = (v_0, ..., v_N), with N+5, N+6 replaced by N+1+σ, N+2+σ ; G[N-k+1] ≥ G_{N-k}
    s = interval(BigFloat, 0)
    for k = 1:N
        s = s + interval(k)*abs(v[k+1])*pown(ν, k)*G[N-k+1]
    end
    return ν*s/interval((N+1)*(N+1+σ)) + abs(μ)*(pown(ν, N+1)*abs(v[N])/interval((N+1)*(N+1+σ)) + pown(ν, N+2)*abs(v[N+1])/interval((N+2)*(N+2+σ)))
end

function solution_operators(φa, ψa, ψb, φb, Λ, ν, N)
    #### (U, rU, V, rV) of Sa = 𝒮^a(1,λ), Ca10 = 𝒞^a(1,0,λ), Ca01 = 𝒞^a(0,1,λ), Cb10, Cb01, Sb = 𝒮^b(1,λ), for λ ∈ Λ
    return solution_operator(φa, Λ, 1, 0, ν, N), solution_operator(ψa, Λ, 1, 0, ν, N), solution_operator(ψa, Λ, 0, 1, ν, N),
           solution_operator(ψb, Λ, 1, 0, ν, N), solution_operator(ψb, Λ, 0, 1, ν, N), solution_operator(φb, Λ, 1, 0, ν, N)
end

function values_at(op, y, h, ν)
    #### (U*, (U*)', ∂λU*, ∂λ(U*)') at y = ±1, ' = d/dx = h^{-1}d/dy, by (A.11) with c_1 = 1/(e ln ν)
    U, rU, V, rV = op
    c1 = inv(exp(interval(BigFloat, 1))*log(ν))
    return (U(y) + rU*interval(-1, 1), (differentiate(U)(y) + c1*rU*interval(-1, 1))/h,
            V(y) + rV*interval(-1, 1), (differentiate(V)(y) + c1*rV*interval(-1, 1))/h)
end



############ Sections A.1 and A.3.1 : the plateau [x2,x3], u_λ, v_λ and their λ-derivatives, with the Bessel functions of Arblib ############

function besselJ(n, z)
    #### J_n(z) at an exact point z > 0 : Arb ball at 128 + 2z bits (the series loses about 1.44z bits), as an Interval{BigFloat}
    p = precision(BigFloat) + 2*ceil(Int, Float64(z))
    return interval(Arblib.hypgeom_bessel_j!(Arb(prec = p), Arb(n, prec = p), Arb(z, prec = p); prec = p))
end

function besselY(n, z)
    #### Y_n(z) at an exact point z > 0, as besselJ
    p = precision(BigFloat) + 2*ceil(Int, Float64(z))
    return interval(Arblib.hypgeom_bessel_y!(Arb(prec = p), Arb(n, prec = p), Arb(z, prec = p); prec = p))
end

#### M_m(z) = √(J_m(z)² + Y_m(z)²), decreasing on (0,∞) (Nicholson's formula) : |J_m(ζ)|, |Y_m(ζ)| ≤ M_m(z) for ζ ≥ z, M_{-m} = M_m
nicholson(m, z) = sqrt(besselJ(abs(m), z)^2 + besselY(abs(m), z)^2)

function bessel_ball(B, n, z)
    #### B_n(ζ) for every ζ in the interval z ⊂ (0,∞), B = besselJ or besselY : mean value theorem from the midpoint (companion (11))
    inf(z) > 0 || error("The argument of the Bessel function is not positive : the proof is not conclusive.")
    return B(n, mid(z)) + (z - interval(mid(z)))*interval(-1, 1)*(nicholson(n-1, inf(z)) + nicholson(n+1, inf(z)))/interval(2)
end

function plateau(radial, x, λ)
    #### (u, u', ∂λu, ∂λu') and (v, v', ∂λv, ∂λv') at x, k = √λ : u = J_1(kx)/x, v = Y_1(kx)/x (ρ = x³P), u = cos kx, v = sin kx (ρ = P)
    k = sqrt(λ) ; z = k*x
    if radial
        J = [bessel_ball(besselJ, n, z) for n = 0:3] ; Y = [bessel_ball(besselY, n, z) for n = 0:3]     #### J[n+1] = J_n(kx)
        u = (J[2]/x, -k*J[3]/x, (J[1] - J[3])/(interval(4)*k), -J[3]/(interval(2)*k*x) - (J[2] - J[4])/interval(4))
        v = (Y[2]/x, -k*Y[3]/x, (Y[1] - Y[3])/(interval(4)*k), -Y[3]/(interval(2)*k*x) - (Y[2] - Y[4])/interval(4))
    else
        s = sin(z) ; c = cos(z)
        u = (c, -k*s, -x*s/(interval(2)*k), -s/(interval(2)*k) - x*c/interval(2))
        v = (s, k*c, x*c/(interval(2)*k), c/(interval(2)*k) - x*s/interval(2))
    end
    return u, v
end

function plateau_bounds(radial, x2, x3, Λ)
    #### S_i ≥ |u^(i)|, |v^(i)| and D_i ≥ |∂λu^(i)|, |∂λv^(i)| on [x2,x3] × Λ, i = 0, 1, 2 (Section A.4, companion (12)) :
    #### u'' = -(3/x)u' - λu (ρ = x³P), u'' = -λu (ρ = P), ∂λu'' = -(3/x)∂λu' - u - λ∂λu ; |J_m(kx)|, |Y_m(kx)| ≤ M_m(k x2)
    klo = interval(inf(sqrt(Λ))) ; khi = interval(sup(sqrt(Λ))) ; λhi = interval(sup(Λ))
    if radial
        x = interval(inf(x2))
        M = [nicholson(m, inf(klo*x)) for m = 0:3]                                #### M[m+1] = M_m(k x2)
        S0 = M[2]/x ; S1 = khi*M[3]/x ; S2 = interval(3)*S1/x + λhi*S0
        D0 = (M[1] + M[3])/(interval(4)*klo) ; D1 = M[3]/(interval(2)*klo*x) + (M[2] + M[4])/interval(4) ; D2 = interval(3)*D1/x + S0 + λhi*D0
    else
        X = interval(max(sup(abs(x2)), sup(abs(x3))))                            #### |x| ≤ X on the plateau
        S0 = interval(1) ; S1 = khi ; S2 = λhi
        D0 = X/(interval(2)*klo) ; D1 = inv(interval(2)*klo) + X/interval(2) ; D2 = S0 + λhi*D0
    end
    return [S0, S1, S2], [D0, D1, D2]
end



############ Section A.2 : L(λ) and ∂λL(λ) of (A.7), g(θ) = L(λ)θ̃ - e_9, Dg(θ) = (∂λL(λ)θ̃ | L(λ)), for all λ ∈ Λ ############

function computation_L(φa, ψa, ψb, φb, ops, radial, x2, x3, Λ, ν)
    #### rows 1-4 : the jumps of φ at x1, x2, x3, x4, rows 5-8 : the jumps of φ', row 9 : a0 - 1 (g_1, ..., g_9) ;
    #### columns θ̃ = (a0, ca0, ca1, A, B, cb0, cb1, b0) ; the second entry of a piece is h
    Sa, Ca10, Ca01, Cb10, Cb01, Sb = ops
    yR = interval(BigFloat, 1)                                                    #### y = 1, and y = -1 for ψa at x1, ψb at x3
    sa = values_at(Sa, yR, φa[2], ν) ; sb = values_at(Sb, yR, φb[2], ν)
    ca10_1 = values_at(Ca10, -yR, ψa[2], ν) ; ca01_1 = values_at(Ca01, -yR, ψa[2], ν)
    ca10_2 = values_at(Ca10, yR, ψa[2], ν)  ; ca01_2 = values_at(Ca01, yR, ψa[2], ν)
    cb10_3 = values_at(Cb10, -yR, ψb[2], ν) ; cb01_3 = values_at(Cb01, -yR, ψb[2], ν)
    cb10_4 = values_at(Cb10, yR, ψb[2], ν)  ; cb01_4 = values_at(Cb01, yR, ψb[2], ν)
    u2, v2 = plateau(radial, x2, Λ) ; u3, v3 = plateau(radial, x3, Λ)
    L = fill(interval(BigFloat, 0), 9, 8) ; dL = fill(interval(BigFloat, 0), 9, 8)
    for k = 0:1                                                                   #### k = 0 : values, rows 1-4 ; k = 1 : derivatives, rows 5-8
        i = 4k ; iL = k+1 ; idL = k+3                                             #### entries iL of the 4-tuples for L, idL for ∂λL
        L[i+1, 1] = sa[iL] ; L[i+1, 2] = -ca10_1[iL] ; L[i+1, 3] = -ca01_1[iL]                                  #### x1 : a0 Sa - ca0 Ca10 - ca1 Ca01
        L[i+2, 2] = ca10_2[iL] ; L[i+2, 3] = ca01_2[iL] ; L[i+2, 4] = -u2[iL] ; L[i+2, 5] = -v2[iL]              #### x2 : ca0 Ca10 + ca1 Ca01 - A u - B v
        L[i+3, 4] = u3[iL] ; L[i+3, 5] = v3[iL] ; L[i+3, 6] = -cb10_3[iL] ; L[i+3, 7] = -cb01_3[iL]              #### x3 : A u + B v - cb0 Cb10 - cb1 Cb01
        L[i+4, 6] = cb10_4[iL] ; L[i+4, 7] = cb01_4[iL] ; L[i+4, 8] = -sb[iL]                                    #### x4 : cb0 Cb10 + cb1 Cb01 - b0 Sb
        dL[i+1, 1] = sa[idL] ; dL[i+1, 2] = -ca10_1[idL] ; dL[i+1, 3] = -ca01_1[idL]
        dL[i+2, 2] = ca10_2[idL] ; dL[i+2, 3] = ca01_2[idL] ; dL[i+2, 4] = -u2[idL] ; dL[i+2, 5] = -v2[idL]
        dL[i+3, 4] = u3[idL] ; dL[i+3, 5] = v3[idL] ; dL[i+3, 6] = -cb10_3[idL] ; dL[i+3, 7] = -cb01_3[idL]
        dL[i+4, 6] = cb10_4[idL] ; dL[i+4, 7] = cb01_4[idL] ; dL[i+4, 8] = -sb[idL]
    end
    L[9, 1] = interval(BigFloat, 1)
    return L, dL
end



############ Section A.4 : enclosure of a polynomial and of its derivative on an interval, number of zeros of φ ############

function enclosure(a, ylo, yhi)
    #### a(y) ∈ Σ_m b_m([ylo,yhi] - y0)^m, b_m = a^(m)(y0)/m! by Horner's scheme at the midpoint y0, then interval Horner (companion (15))
    N = length(coefficients(a)) - 1
    y0 = (ylo + yhi)/2
    b = copy(coefficients(a))
    for m = 0:N-1
        for n = N-1:-1:m
            b[n+1] = b[n+1] + interval(y0)*b[n+2]
        end
    end
    s = interval(ylo, yhi) - interval(y0)
    value = b[N+1] ; deriv = interval(N)*b[N+1]
    for m = N-1:-1:1
        value = value*s + b[m+1] ; deriv = deriv*s + interval(m)*b[m+1]
    end
    return value*s + b[1], deriv
end

function series_on(piece, p, X)
    #### (φ̃(x), φ̃'(x)) for x in the interval X on a series piece, φ̃ = the polynomial p in y = (x - c)/h, ' = d/dx
    c, h = piece
    Y = (X - c)/h
    value, deriv = enclosure(p, inf(Y), sup(Y))
    return value, deriv/h
end

function plateau_on(radial, λ, A, B, X)
    #### (φ̃(x), φ̃'(x)) for x in the interval X on the plateau, φ̃ = A u_λ + B v_λ
    u, v = plateau(radial, X, λ)
    return A*u[1] + B*v[1], A*u[2] + B*v[2]
end

function classify(f, x1, x2, ε, depth)
    #### [x1,x2] is halved until f = (φ̃, φ̃') shows φ̃ > ε (+), φ̃ < -ε (-) or |φ̃'| > ε (0, the sign of φ̃' recorded) ; also sup|φ̃|
    value, deriv = f(interval(x1, x2))
    (isguaranteed(value) && isguaranteed(deriv)) || error("An enclosure of the zero count is not guaranteed : the proof is not conclusive.")
    if inf(value) > ε
        return [('+', 0, x1, x2)], sup(abs(value))
    elseif sup(value) < -ε
        return [('-', 0, x1, x2)], sup(abs(value))
    elseif inf(deriv) > ε
        return [('0', 1, x1, x2)], sup(abs(value))
    elseif sup(deriv) < -ε
        return [('0', -1, x1, x2)], sup(abs(value))
    elseif depth < 40
        m = (x1 + x2)/2
        l1, s1 = classify(f, x1, m, ε, depth+1) ; l2, s2 = classify(f, m, x2, ε, depth+1)
        return [l1; l2], max(s1, s2)
    else
        error("The zero count is undecided near x = $(Float64(x1)) : the proof is not conclusive.")
    end
end

function count_zeros(nodes, radial, pieces, φ̃, λ̄, Ā, B̄, ε)
    #### every f ∈ C¹([a,b]) with ‖f - φ̃‖_{C¹} ≤ ε has the sign of φ̃ on the intervals (+), (-), and is strictly monotone on each maximal
    #### run of intervals (0) : one zero in each run between a (+) and a (-) interval, no other zero (Section A.4) ; about √λ̄ intervals
    #### per unit length at the start. Returns the number of zeros, the intervals containing them, and sup|φ̃| on the four series pieces.
    a, x1, x2, x3, x4, b = nodes
    k̄ = sqrt(Float64(mid(λ̄)))
    leaves = [] ; supφ̃ = zeros(BigFloat, 4)
    sweep = [(a, x1, X -> series_on(pieces[1], φ̃[1], X), 1), (x1, x2, X -> series_on(pieces[2], φ̃[2], X), 2),
             (x2, x3, X -> plateau_on(radial, λ̄, Ā, B̄, X), 0),
             (x3, x4, X -> series_on(pieces[3], φ̃[3], X), 3), (x4, b, X -> series_on(pieces[4], φ̃[4], X), 4)]
    for (lo, hi, f, j) in sweep
        n = max(2, ceil(Int, Float64(sup(hi - lo))*k̄))
        grid = [inf(lo) + (sup(hi) - inf(lo))*i/n for i = 0:n] ; grid[1] = inf(lo) ; grid[n+1] = sup(hi)
        for i = 1:n
            l, s = classify(f, grid[i], grid[i+1], ε, 0)
            append!(leaves, l)
            if j > 0
                supφ̃[j] = max(supφ̃[j], s)
            end
        end
    end
    runs = []                                                                     #### consecutive intervals (0) are merged
    for (kind, sgn, y1, y2) in leaves
        if kind == '0' && !isempty(runs) && runs[end][1] == '0'
            runs[end][2] == sgn || error("φ̃' changes sign inside a run of the zero count : the proof is not conclusive.")
            runs[end] = ('0', sgn, runs[end][3], y2)
        else
            push!(runs, (kind, sgn, y1, y2))
        end
    end
    (runs[1][1] != '0' && runs[end][1] != '0') || error("An interval (0) of the zero count touches a or b : the proof is not conclusive.")
    for i = 2:length(runs)
        (runs[i][1] == '0' || runs[i-1][1] == '0' || runs[i][1] == runs[i-1][1]) ||
            error("φ̃ changes sign without a run of the zero count : the proof is not conclusive.")
    end
    zeros_in = [(runs[i][3], runs[i][4]) for i = 2:length(runs)-1 if runs[i][1] == '0' && runs[i-1][1] != runs[i+1][1]]
    return length(zeros_in), zeros_in, supφ̃
end



############ Section A.4 : normalization 𝒩 = ∫_a^b ρφ² dx and the part of parity (-1)^m of ψ̃ ############

function computation_𝒩(pieces, φ̃, e, supφ̃, radial, x2, x3, Λ, A, B)
    #### series piece : |h|∫Rφ̃² dy exactly (polynomial), plus |∫R(φ² - φ̃²)| ≤ (2e sup|φ̃| + e²)∫R since ρ ≥ 0 and |φ - φ̃| ≤ e ;
    #### plateau (k = √λ, λ ∈ Λ, A, B ∋ θ*_5, θ*_6) : ∫x Z_1(kx)² dx = x²/2 (Z_1² - Z_0 Z_2)(kx), Z_n = A J_n + B Y_n (Lommel, ρ = x³P),
    #### ∫(A cos kx + B sin kx)² dx = (A² + B²)x/2 + (A² - B²)sin(2kx)/(4k) - AB cos(2kx)/(2k) (ρ = P)
    𝒩 = interval(BigFloat, 0)
    for j = 1:4
        c, h, σ, R = pieces[j]
        P = R*φ̃[j]*φ̃[j]
        dR = length(coefficients(R)) - 1 ; dP = length(coefficients(P)) - 1
        if σ == 4                                                                 #### φa, φb : y ∈ [0,1]
            intR = sum(R[k]/interval(k+1) for k = 0:dR) ; intP = sum(P[k]/interval(k+1) for k = 0:dP)
        else                                                                      #### ψa, ψb : y ∈ [-1,1], odd powers integrate to 0
            intR = sum(interval(2)*R[k]/interval(k+1) for k = 0:2:dR) ; intP = sum(interval(2)*P[k]/interval(k+1) for k = 0:2:dP)
        end
        𝒩 = 𝒩 + abs(h)*(intP + interval(-1, 1)*(interval(2)*e[j]*interval(supφ̃[j]) + e[j]^2)*intR)
    end
    k = sqrt(Λ)
    function primitive(x)
        #### a primitive of ρ(A u_λ + B v_λ)² on the plateau
        if radial
            Z0 = A*bessel_ball(besselJ, 0, k*x) + B*bessel_ball(besselY, 0, k*x)
            Z1 = A*bessel_ball(besselJ, 1, k*x) + B*bessel_ball(besselY, 1, k*x)
            Z2 = A*bessel_ball(besselJ, 2, k*x) + B*bessel_ball(besselY, 2, k*x)
            return x^2/interval(2)*(Z1^2 - Z0*Z2)
        else
            return (A^2 + B^2)*x/interval(2) + (A^2 - B^2)*sin(interval(2)*k*x)/(interval(4)*k) - A*B*cos(interval(2)*k*x)/(interval(2)*k)
        end
    end
    return 𝒩 + primitive(x3) - primitive(x2)
end

function symmetrise(ψ̃, Ã, B̃, parity)
    #### the part of parity (-1)^m of ψ̃ (even weight, Lemma 4.13 (v)) : the image of φa is φb with the same y, the one of ψa is ψb with
    #### y ↦ -y ; on the plateau (ρ = P) the part of Ã cos kx + B̃ sin kx is Ã cos kx (parity 1) or B̃ sin kx (parity -1)
    N = length(coefficients(ψ̃[1])) - 1
    flip(p) = Sequence(Taylor(N), [interval((-1)^n)*coefficients(p)[n+1] for n = 0:N])
    φa = (ψ̃[1] + interval(parity)*ψ̃[4])/interval(2) ; ψa = (ψ̃[2] + interval(parity)*flip(ψ̃[3]))/interval(2)
    return (φa, ψa, interval(parity)*flip(ψa), interval(parity)*φa), (parity == 1 ? Ã : interval(BigFloat, 0)), (parity == 1 ? interval(BigFloat, 0) : B̃)
end
