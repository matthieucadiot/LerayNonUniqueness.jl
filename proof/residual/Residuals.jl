# Residuals.jl : the residual numerators on a cell (Subsections "Polynomial numerators in the interior" and
# "Computation on a bounded domain").
#
# On a cell, a field is a named tuple (f, u, v, uu, vv) of models of the field and of its derivatives
# ∂_u, ∂_v, ∂_u^2, ∂_v^2. The coefficients of the operators (u:Luv), (u:K), (u:Lpuv), multiplied by D = A^7 B^7,
# are the polynomials of A = 1+u and B = 1+v collected in cell_polynomials:
#   D L_uv = B^7 [3A(A^4-1) ∂_u^2 + 3(5A^4+3) ∂_u] + A^7 [3B(B^4-1) ∂_v^2 + 9(B^4+1) ∂_v]        (for L̃_uv: 3(3-B^4) ∂_v)
#   D K(u) = A^4(A^4-1)/4 B^7,   D K(v) = A^7 B^4(B^4-1)/4.

# the polynomials A^7, 3A(A^4-1), 3(5A^4+3), A^4(A^4-1)/4 of A = 1+u on the u-cell [a,b] (thin intervals), as exact models of x
function cell_polynomials_u(a, b, N)
    A = zeros(Interval{BigFloat}, Taylor(1) ⊗ Taylor(0))
    A[(0,0)] = I"1" + (a + b)/I"2"
    A[(1,0)] = (b - a)/I"2"
    A = Model(A, N)
    A4 = A^4
    return (A = A, D = A^7, L2 = 3*A*(A4 - I"1"), L1 = 3*(5*A4 + I"3"), K = I"0.25"*A4*(A4 - I"1"))
end

# the polynomials B^7, 3B(B^4-1), 9(B^4+1), 3(3-B^4), B^4(B^4-1)/4 of B = 1+v on the v-cell [c,d], as exact models of y
function cell_polynomials_v(c, d, N)
    B = zeros(Interval{BigFloat}, Taylor(0) ⊗ Taylor(1))
    B[(0,0)] = I"1" + (c + d)/I"2"
    B[(0,1)] = (d - c)/I"2"
    B = Model(B, N)
    B4 = B^4
    return (B = B, D = B^7, L2 = 3*B*(B4 - I"1"), L1 = 9*(B4 + I"1"), L̃1 = 3*(I"3" - B4), K = I"0.25"*B4*(B4 - I"1"))
end

# N_Ψ = D (L_uv Ψ - Ω), first equation of (u:nl)
function numerator_Ψ(Ψ, Ω, a, b)
    return b.D*(a.L2*Ψ.uu + a.L1*Ψ.u) + a.D*(b.L2*Ψ.vv + b.L1*Ψ.v) - a.D*b.D*Ω.f
end

# N_Ω = D (L_uv Ω + [K(u)∂_u + K(v)∂_v + 2] Ω + 2 𝒩[Ψ,Ω]), second equation of (u:nl), with 𝒩 of (u:N)
function numerator_Ω(Ψ, Ω, a, b)
    LΩ = b.D*(a.L2*Ω.uu + a.L1*Ω.u) + a.D*(b.L2*Ω.vv + b.L1*Ω.v)
    𝒩 = 2*a.K*b.K*(Ψ.u*Ω.v - Ψ.v*Ω.u) - a.K*b.D*(Ψ.f*Ω.u - Ψ.u*Ω.f) + 2*a.D*b.K*Ψ.f*Ω.v + a.D*b.D*Ψ.f*Ω.f
    return LΩ + a.K*b.D*Ω.u + a.D*b.K*Ω.v + 2*a.D*b.D*Ω.f + 2*𝒩
end

# Ñ_Ψ = D (L̃_uv Ψ̃ - Ω̃), first equation of (u:eig)
function numerator_Ψ̃(Ψ̃, Ω̃, a, b)
    return b.D*(a.L2*Ψ̃.uu + a.L1*Ψ̃.u) + a.D*(b.L2*Ψ̃.vv + b.L̃1*Ψ̃.v) - a.D*b.D*Ω̃.f
end

# Ñ_Ω = D (L̃_uv Ω̃ + [K(u)∂_u + K(v)∂_v + 3/2 - λ̃] Ω̃ + 2 𝒩̃_1(Ψ,Ω̃) + 2 𝒩̃_2(Ψ̃,Ω)), second equation of (u:eig), (u:halves)
function numerator_Ω̃(Ψ, Ω, Ψ̃, Ω̃, λ̃, a, b)
    LΩ̃ = b.D*(a.L2*Ω̃.uu + a.L1*Ω̃.u) + a.D*(b.L2*Ω̃.vv + b.L̃1*Ω̃.v)
    𝒩1 = 2*a.K*b.K*(Ψ.u*Ω̃.v - Ψ.v*Ω̃.u) - a.K*b.D*Ψ.f*Ω̃.u + 2*a.D*b.K*Ψ.f*Ω̃.v
    𝒩2 = 2*a.K*b.K*(Ψ̃.u*Ω.v - Ψ̃.v*Ω.u) + a.K*b.D*Ψ̃.u*Ω.f + 2*a.D*b.K*Ψ̃.f*Ω.v + a.D*b.D*Ψ̃.f*Ω.f
    return LΩ̃ + a.K*b.D*Ω̃.u + a.D*b.K*Ω̃.v + (I"1.5" - λ̃)*a.D*b.D*Ω̃.f + 2*𝒩1 + 2*𝒩2
end

# ---------------------------------------------------------------------------------------------------------------
# The gluing region S_band = S_L \ S_{L-d}

# χ_0(u) = 1 - ϑ((u-L+d)/d) of (res:cutoff), with ϑ of (sec3eq5), on the u-cell [a,b] (thin intervals) as an exact model
# of x (and χ_0 = 1 if the cell is contained in [0,L-d])
function cutoff_u(a, b, L, d, N)
    if sup(b) ≤ L - d
        return Model(Sequence(Taylor(0) ⊗ Taylor(0), [big(I"1")]), N)
    end
    inf(a) ≥ L - d && sup(b) ≤ L || error("the cell [$a,$b] is not contained in [L-d,L]")
    τ = zeros(Interval{BigFloat}, Taylor(1) ⊗ Taylor(0))
    τ[(0,0)] = ((a + b)/I"2" - interval(L - d))/interval(d)
    τ[(1,0)] = (b - a)/(I"2"*interval(d))
    τ = Model(τ, N)
    ϑ = 126*τ^5 - 420*τ^6 + 540*τ^7 - 315*τ^8 + 70*τ^9
    return I"1" - ϑ
end

# the same in the variable y of the v-cell [c,e]
function cutoff_v(c, e, L, d, N)
    if sup(e) ≤ L - d
        return Model(Sequence(Taylor(0) ⊗ Taylor(0), [big(I"1")]), N)
    end
    inf(c) ≥ L - d && sup(e) ≤ L || error("the cell [$c,$e] is not contained in [L-d,L]")
    τ = zeros(Interval{BigFloat}, Taylor(0) ⊗ Taylor(1))
    τ[(0,0)] = ((c + e)/I"2" - interval(L - d))/interval(d)
    τ[(0,1)] = (e - c)/(I"2"*interval(d))
    τ = Model(τ, N)
    ϑ = 126*τ^5 - 420*τ^6 + 540*τ^7 - 315*τ^8 + 70*τ^9
    return I"1" - ϑ
end

# χ = χ_0(u) χ_0(v) and its derivatives on the cell, as exact models, from χ_0(u), χ_0(v) and the half-widths h_u, h_v
function cutoff(χu, χv, hu, hv)
    inf(χu.ϱ) == sup(χu.ϱ) == 0 && inf(χv.ϱ) == sup(χv.ϱ) == 0 || error("the cut-off polynomials must be exact before differentiation")
    inf(hu) > 0 && inf(hv) > 0 || error("the cell half-widths must be positive")
    χu_x = Model(differentiate(χu.Π, (1,0))/hu, χu.N)
    χu_xx = Model(differentiate(χu_x.Π, (1,0))/hu, χu.N)
    χv_y = Model(differentiate(χv.Π, (0,1))/hv, χv.N)
    χv_yy = Model(differentiate(χv_y.Π, (0,1))/hv, χv.N)
    return (f = χu*χv, u = χu_x*χv, v = χu*χv_y, uu = χu_xx*χv, vv = χu*χv_yy)
end

# the glued field f = f_ext + χ (f_in - f_ext) and its derivatives by the product rule, (X:glue), (X:glueder)
function glued_field(f_in, f_ext, χ)
    g = f_in.f - f_ext.f
    g_u = f_in.u - f_ext.u
    g_v = f_in.v - f_ext.v
    g_uu = f_in.uu - f_ext.uu
    g_vv = f_in.vv - f_ext.vv
    return (f = f_ext.f + χ.f*g,
            u = f_ext.u + χ.u*g + χ.f*g_u,
            v = f_ext.v + χ.v*g + χ.f*g_v,
            uu = f_ext.uu + χ.uu*g + 2*χ.u*g_u + χ.f*g_uu,
            vv = f_ext.vv + χ.vv*g + 2*χ.v*g_v + χ.f*g_vv)
end

# the models of x = r^2/R^2, t = 1/R and ℓ = log R of (cap:asymptotic-polynomial) on the cell, from
# R^2 = U_0(u) + U_0(v), t = (R^2)^{-1/2}, x = U_0(u) t^2, ℓ = (log R^2)/2, with the models A, B of 1+u, 1+v of the
# cell polynomials (truncated at the bidegree N of the gluing region)
function cell_variables(A, B, N)
    A = Model(A.Π, A.ϱ, N)
    B = Model(B.Π, B.ϱ, N)
    U0u = (A^4 - I"1")/I"12"
    U0v = (B^4 - I"1")/I"12"
    R2 = U0u + U0v
    T = R2^I"-0.5"
    X = U0u*(T*T)
    L = I"0.5"*log(R2)
    return X, T, L
end

# the models of the monomials x^k t^h ℓ^j, k ≤ K, h ≤ H, j ≤ J, on the cell
function monomials(X, T, L, K, H, J)
    XTL = Array{Model}(undef, K+1, H+1, J+1)
    XTL[1,1,1] = Model(one(X.Π), X.N)
    for k = 1:K
        XTL[k+1,1,1] = XTL[k,1,1]*X
    end
    for k = 0:K, h = 1:H
        XTL[k+1,h+1,1] = XTL[k+1,h,1]*T
    end
    for k = 0:K, h = 0:H, j = 1:J
        XTL[k+1,h+1,j+1] = XTL[k+1,h+1,j]*L
    end
    return XTL
end

# the model of F(x,t,ℓ) = Σ F_{khj} x^k t^h ℓ^j on the cell, for a far-field polynomial F
function far_evaluate(F, XTL)
    G = I"0"*XTL[1,1,1]
    for (k,h,j) in indices(space(F))
        c = F[(k,h,j)]
        if !(inf(c) == 0 && sup(c) == 0)
            G = G + c*XTL[k+1,h+1,j+1]
        end
    end
    return G
end

# the asymptotic field on a band cell with its (u,v)-derivatives, (X:uderiv):
#   ∂_u f = (A^3/6) (∂_r f)/r,   ∂_u^2 f = (A^2/2) (∂_r f)/r + (A^6/36) (∂_r/r)^2 f,   and the same in v with B,
# where the far-field polynomials of (∂_r f)/r = t^2 𝒜_r f, (∂_r/r)^2 f = t^4(𝒜_r^2 f - 2 𝒜_r f), and the z-analogues,
# are the entries r1, r2, z1, z2 of F (Asymptotics.jl, far_field_derivatives)
function far_field_on_cell(F, XTL, A, B)
    f = far_evaluate(F.f, XTL)
    fr1 = far_evaluate(F.r1, XTL)
    fr2 = far_evaluate(F.r2, XTL)
    fz1 = far_evaluate(F.z1, XTL)
    fz2 = far_evaluate(F.z2, XTL)
    A2 = A*A; A3 = A2*A; A6 = A3*A3
    B2 = B*B; B3 = B2*B; B6 = B3*B3
    return (f = f,
            u = (one(A.Π[(0,0)])/I"6")*A3*fr1,
            v = (one(B.Π[(0,0)])/I"6")*B3*fz1,
            uu = I"0.5"*A2*fr1 + (one(A.Π[(0,0)])/I"36")*A6*fr2,
            vv = I"0.5"*B2*fz1 + (one(B.Π[(0,0)])/I"36")*B6*fz2)
end
