# Norms.jl : the cell integrals, the coefficient bounds and the assembly of Table certi1 and of (res:main-bound),
# (res:even-total) (Section "Evaluation of the residual norms").
#
# On a cell C = [u_i,u_{i+1}] × [v_j,v_{j+1}] the measure (X:measure) and the factors z^2 (odd fields), r^2 give
#   z^2 J du dv = (π²/9) U_0(u) A^3 B^3 U_0(v)^{1/2} du dv,      J du dv = (π²/9) U_0(u) A^3 B^3 U_0(v)^{-1/2} du dv,
# and a residual E = A^{-7} B^{-7} N brings the factor A^{-14} B^{-14}. All the densities are products ρ(x) μ(y) of
# one-dimensional models (prepared by density()), and the cell integral of N^2 ρ μ is computed by cell_integral.

U0(x) = ((x + I"1")^I"4" - I"1")/I"12"                 # for an interval x

# everything that depends only on the u-cell [U_i,U_{i+1}]: the local basis, the cell polynomials, the cut-off, the
# u-densities (π²/9) U_0 A^3, (π²/9) U_0^2 A^3, (π²/9) U_0 A^{-11}, (π²/9) U_0^2 A^{-11}, and the bounds on r = U_0(u)^{1/2}
function column_data(i, data, par)
    a = data.U[i]; b = data.U[i+1]                                     # thin 128-bit intervals
    af = a; bf = b                                                   # retain the 128-bit endpoints
    A = zeros(Interval{Float64}, Taylor(1) ⊗ Taylor(0))
    A[(0,0)] = I"1" + (af + bf)/I"2"
    A[(1,0)] = (bf - af)/I"2"
    A = Model(A, par.N_ρ)
    U0u = (A^4 - I"1")/I"12"
    ρJ = RESPONSE_UNIT_PI64*RESPONSE_UNIT_PI64/I"9"*U0u*A^3
    ρres = RESPONSE_UNIT_PI64*RESPONSE_UNIT_PI64/I"9"*U0u*A^I"-11"
    return (a = a, b = b, h = interval(Float64,(b - a)/I"2"), hf = interval(Float64,(bf - af)/I"2"),
            k = knot_interval(data.tu, a), B = local_basis(data.tu, data.p, knot_interval(data.tu, a), a, b),
            poly = cell_polynomials_u(a, b, par.N_in), χ = cutoff_u(a, b, data.L, data.d, par.N_in), interior = sup(b) ≤ data.L - data.d,
            ρJ = density(ρJ, par.dmom), ρJr = density(U0u*ρJ, par.dmom), ρres = density(ρres, par.dmom), ρresr = density(U0u*ρres, par.dmom),
            U0a = interval(Float64,U0(af)), U0b = interval(Float64,U0(bf)), r = interval(Float64,sqrt(U0(bf))))
end

# everything that depends only on the v-cell [V_j,V_{j+1}]: as above, with the v-densities B^3 U_0(v)^{±1/2} and
# B^{-11} U_0(v)^{±1/2}, the model of z = U_0(v)^{1/2} for the suprema, and the Jacobian of the variable of integration.
# On the cell touching the axis, V_j = 0, the chart (cap:axis-chart) v = s^2, s = (√e/2)(1+σ), is used: the densities
# (cap:axis-density) are 2 s^2 B^3 H^{1/2}, 2 B^3 H^{-1/2}, 2 s^2 B^{-11} H^{1/2}, 2 B^{-11} H^{-1/2} with B = 1+s^2,
# H = (4+6s^2+4s^4+s^6)/12, the model of z is s H^{1/2}, and the Jacobian is √e/2.
function row_data(j, data, par)
    c = data.V[j]; e = data.V[j+1]                                     # thin 128-bit intervals
    cf = c; ef = e                                                   # retain the 128-bit endpoints
    axis = (sup(c) == 0)
    if axis
        s = zeros(Interval{Float64}, Taylor(0) ⊗ Taylor(1))
        s[(0,0)] = sqrt(ef)/I"2"
        s[(0,1)] = sqrt(ef)/I"2"
        s = Model(s, par.N_ρ)
        s2 = s*s
        B = I"1" + s2
        H = (I"4" + I"6"*s2 + I"4"*s2*s2 + s2*s2*s2)/I"12"
        μJodd = I"2"*s2*B^3*H^I"0.5"
        μJeven = I"2"*B^3*H^I"-0.5"
        μres_odd = I"2"*s2*B^I"-11"*H^I"0.5"
        μres_even = I"2"*B^I"-11"*H^I"-0.5"
        z = s*H^I"0.5"
        jac = interval(Float64,sqrt(ef)/I"2")
        dmom = 2*par.dmom
    else
        B = zeros(Interval{Float64}, Taylor(0) ⊗ Taylor(1))
        B[(0,0)] = I"1" + (cf + ef)/I"2"
        B[(0,1)] = (ef - cf)/I"2"
        B = Model(B, par.N_ρ)
        U0v = (B^4 - I"1")/I"12"
        μJodd = B^3*U0v^I"0.5"
        μJeven = B^3*U0v^I"-0.5"
        μres_odd = B^I"-11"*U0v^I"0.5"
        μres_even = B^I"-11"*U0v^I"-0.5"
        z = U0v^I"0.5"
        jac = interval(Float64,(ef - cf)/I"2")
        dmom = par.dmom
    end
    return (c = c, e = e, h = interval(Float64,(e - c)/I"2"), jac = jac, axis = axis,
            k = knot_interval(data.tv, c), B = local_basis(data.tv, data.p, knot_interval(data.tv, c), c, e),
            poly = cell_polynomials_v(c, e, par.N_in), χ = cutoff_v(c, e, data.L, data.d, par.N_in), interior = sup(e) ≤ data.L - data.d,
            μJodd = density(μJodd, dmom), μJeven = density(μJeven, dmom), μres_odd = density(μres_odd, dmom), μres_even = density(μres_even, dmom),
            z = Model(z.Π, z.ϱ, par.N_in), U0c = interval(Float64,U0(cf)), U0e = interval(Float64,U0(ef)), zsup = interval(Float64,sqrt(U0(ef))))
end

# Φ(R)^2 of (res:Phi), R < R_1
function Φ2(R, R1)
    R = interval(Float64, R); R1 = interval(Float64, R1)
    0 ≤ inf(R) && sup(R) < inf(R1) || error("Φ requires 0 ≤ R < R1")
    S = I"0"
    for j = 0:4
        S = S + interval(Float64, binomial(4,j))*R^interval(4-j)/(interval(Float64, 7-j)*(R1 - R)^interval(7-j))
    end
    return I"3"/(I"8"*RESPONSE_UNIT_PI64*RESPONSE_UNIT_PI64)*S
end

# The quantities of the cell (i,j) (Subsections "Computation on a bounded domain" and "The coefficient bounds"):
#   [1] w̄_C ∫_C E_Ω^2 z^2 J        [2] w̄_C ∫_C E_Ψ^2 z^2 J (cell in S_{84})      [3] ∫_C E_Ψ^2 z^2 J (cell in S_band)
#   [4] w̄_C ∫_C E_Ψ^2 z^2 J (cell in S_band)        [5] ∫_C r^2 E_Ω̃^2 J        [6] ∫_C r^2 E_Ψ̃^2 J
#   [7] ∫_C Ω̃^2 J        [8] ∫_C r^2 Ω̃^2 J        [9], [10] the cell terms of I_o(Ω*,R_0)^2 and I_o(Ω̃*,R_0)^2, (cap:coefficient-integrals)
#   [11] ‖rΩ*‖_{L^∞(C)}        [12] ‖rΩ̃*‖_{L^∞(C)}, (cap:cell-supremum)        [13], [14] the same when C meets {R ≥ R_0},
# where w̄_C = w(r_+,z_+) is the value of the Hardy weight (res:w) at the upper corner of the cell.
function cell_bounds(i, j, data, far, col, row, par)
    precision(BigFloat) == par.input_prec || error("cell ($i,$j): the working precision is not $(par.prec) bits")
    C = col[i]; R = row[j]
    interior = C.interior && R.interior
    N = interior ? par.N_in : par.N_band
    Ψ = field(piece(data.Ψ, C.B, C.k, R.B, R.k, data.p), C.h, R.h, N)
    Ω = field(piece(data.Ω, C.B, C.k, R.B, R.k, data.p), C.h, R.h, N)
    Ψ̃ = field(piece(data.Ψ̃, C.B, C.k, R.B, R.k, data.p), C.h, R.h, N)
    Ω̃ = field(piece(data.Ω̃, C.B, C.k, R.B, R.k, data.p), C.h, R.h, N)
    if !interior                                                   # the glued fields (X:glue)
        X, T, L = cell_variables(C.poly.A, R.poly.B, par.N_band)
        XTL = monomials(X, T, L, far.K, far.H, far.J)
        χ = cutoff(C.χ, R.χ, C.h, R.h)
        Ψ = glued_field(Ψ, far_field_on_cell(far.P, XTL, C.poly.A, R.poly.B), χ)
        Ω = glued_field(Ω, far_field_on_cell(far.O, XTL, C.poly.A, R.poly.B), χ)
        Ψ̃ = glued_field(Ψ̃, far_field_on_cell(far.P̃, XTL, C.poly.A, R.poly.B), χ)
        Ω̃ = glued_field(Ω̃, far_field_on_cell(far.Õ, XTL, C.poly.A, R.poly.B), χ)
    end
    a = C.poly; b = R.poly
    NΨ = numerator_Ψ(Ψ, Ω, a, b)
    NΩ = numerator_Ω(Ψ, Ω, a, b)
    NΨ̃ = numerator_Ψ̃(Ψ̃, Ω̃, a, b)
    NΩ̃ = numerator_Ω̃(Ψ, Ω, Ψ̃, Ω̃, data.λ̃, a, b)
    Ωf = Ω.f
    Ω̃f = Ω̃.f
    if R.axis                                                      # the chart (cap:axis-chart), in the working precision
        NΨ = substitute_axis(NΨ); NΩ = substitute_axis(NΩ); NΨ̃ = substitute_axis(NΨ̃); NΩ̃ = substitute_axis(NΩ̃)
        Ωf = substitute_axis(Ωf); Ω̃f = substitute_axis(Ω̃f)
    end
    # Squares, densities and moments use native Float64 intervals.
    jac = C.hf*R.jac
    # the Hardy weight (res:w) at the upper corner, with k(z)^2 = z(1+z)^2(2+z)/6 of (res:k)
    r2 = C.U0b
    z = R.zsup
    k2 = z*(I"1" + z)*(I"1" + z)*(I"2" + z)/I"6"
    w̄ = r2*k2/(r2 + k2)
    # the cell integrals
    IΩ = jac*cell_integral(NΩ, C.ρres, R.μres_odd)
    IΨ = jac*cell_integral(NΨ, C.ρres, R.μres_odd)
    IΩ̃ = jac*cell_integral(NΩ̃, C.ρresr, R.μres_even)
    IΨ̃ = jac*cell_integral(NΨ̃, C.ρresr, R.μres_even)
    nΩ̃ = jac*cell_integral(Ω̃f, C.ρJ, R.μJeven)
    rΩ̃ = jac*cell_integral(Ω̃f, C.ρJr, R.μJeven)
    nΩ = jac*cell_integral(Ωf, C.ρJ, R.μJodd)
    # the coefficient integrals (cap:coefficient-integrals) on the cells meeting B_{R_0}, R_1 = r_* = U_0(L-d)^{1/2}
    Rm = sqrt(C.U0a + R.U0c)
    Rp = sqrt(C.U0b + R.U0e)
    IoΩ = I"0"; IoΩ̃ = I"0"
    if inf(Rm) < par.R0
        par.R0 < inf(data.r★) || error("R_0 must be smaller than r_* = U_0(L-d)^{1/2}, (res:Phi)")
        Φ² = Φ2(min(interval(par.R0), Rp), data.r★)
        IoΩ = Φ²*(C.r + I"2"*z)*(C.r + I"2"*z)*nΩ
        IoΩ̃ = Φ²*(C.r + I"2"*z)*(C.r + I"2"*z)*nΩ̃
    end
    # the suprema (cap:cell-supremum)
    sΩ = C.r*sup_bound(R.z*Ωf)
    sΩ̃ = C.r*sup_bound(Ω̃f)
    distant = sup(Rp) ≥ par.R0
    o = I"0"
    return [w̄*IΩ, interior ? w̄*IΨ : o, interior ? o : IΨ, interior ? o : w̄*IΨ, IΩ̃, IΨ̃, nΩ̃, rΩ̃, IoΩ, IoΩ̃, sΩ, sΩ̃, distant ? sΩ : o, distant ? sΩ̃ : o]
end

# The bounds of Table certi1 and of Lemma res:bounds-lemma from the cell quantities S (14 × number of columns) and
# the exterior contributions of Subsection "Computation on an exterior domain", (cap:exterior-norms), with L_e = U_0(L)^{1/2}
function assemble(S, far, data, par)
    size(S) == (14, length(data.U)-1) || error("the cell array has the wrong dimensions")
    all(x -> isguaranteed(x) && isfinite(inf(x)) && isfinite(sup(x)) && inf(x) ≥ 0, S) || error("the cell quantities are not nonnegative guaranteed finite intervals")
    L = sqrt(U0(interval(Float64, data.L)))
    Σ = [sum(S[k,:]) for k = 1:10]
    EΩ = far.EΩ; EΩ̃ = far.EΩ̃; O = far.O.f; Õ = far.Õ.f
    wRΩ = sqrt(Σ[1] + exterior_integral(EΩ, 2, 2, L))             # ‖√w ℛ_Ω‖, w ≤ r^2 and E_Ω = z (E_Ω/z) in S_ext
    wRΨin = sqrt(Σ[2])                                             # ‖√w ℛ_Ψ^in‖
    RΨband = sqrt(Σ[3])                                            # ‖ℛ_Ψ^band‖
    wRΨband = sqrt(Σ[4])                                           # ‖√w ℛ_Ψ^band‖
    rRΩ̃ = sqrt(Σ[5] + exterior_integral(EΩ̃, 2, 0, L))             # ‖r ℛ_Ω̃‖
    rRΨ̃ = sqrt(Σ[6])                                               # ‖r ℛ_Ψ̃‖
    eN = exterior_integral(Õ, 0, 0, L); eR = exterior_integral(Õ, 2, 0, L)   # on {R ≥ L_e} ⊇ S_ext: upper bounds only
    nΩ̃ = interval(inf(sqrt(Σ[7])), sup(sqrt(Σ[7] + eN)))           # ‖Ω̃*‖, the lower endpoint from the cells alone
    rΩ̃ = interval(inf(sqrt(Σ[8])), sup(sqrt(Σ[8] + eR)))           # ‖r Ω̃*‖
    IoΩ = sqrt(Σ[9]); IoΩ̃ = sqrt(Σ[10])                            # I_o(Ω*,R_0), I_o(Ω̃*,R_0)
    sΩ = max(maximum(S[11,:]), exterior_sup(O, 1, L))              # ‖rΩ*‖_∞
    sΩ̃ = max(maximum(S[12,:]), exterior_sup(Õ, 0, L))              # ‖rΩ̃*‖_∞
    sΩfar = max(maximum(S[13,:]), exterior_sup(O, 1, L))           # ‖rΩ*‖_{L^∞(ᶜB_{R_0})}
    sΩ̃far = max(maximum(S[14,:]), exterior_sup(Õ, 0, L))           # ‖rΩ̃*‖_{L^∞(ᶜB_{R_0})}
    # (res:main-bound) and (res:even-total)
    ℛ = wRΩ + sΩ*wRΨin + sqrt(IoΩ*IoΩ*RΨband*RΨband + sΩfar*sΩfar*wRΨband*wRΨband)
    ℛ̃ = rRΩ̃ + sΩ*rRΨ̃ + sΩ̃*wRΨin + sqrt(IoΩ̃*IoΩ̃*RΨband*RΨband + sΩ̃far*sΩ̃far*wRΨband*wRΨband)
    display("‖√w ℛ_Ω‖, ‖√w ℛ_Ψ^in‖, ‖ℛ_Ψ^band‖, ‖√w ℛ_Ψ^band‖")
    display([wRΩ, wRΨin, RΨband, wRΨband])
    display("‖r ℛ_Ω̃‖, ‖r ℛ_Ψ̃‖")
    display([rRΩ̃, rRΨ̃])
    display("‖Ω̃*‖_{L^2}, ‖r Ω̃*‖_{L^2}")
    display([nΩ̃, rΩ̃])
    display("‖rΩ*‖_∞, ‖rΩ*‖_{L^∞(ᶜB_{R_0})}, I_o(Ω*,R_0)")
    display([sΩ, sΩfar, IoΩ])
    display("‖rΩ̃*‖_∞, ‖rΩ̃*‖_{L^∞(ᶜB_{R_0})}, I_o(Ω̃*,R_0)")
    display([sΩ̃, sΩ̃far, IoΩ̃])
    display("‖ℛ‖_{Z'} ≤ ")
    display(sup(ℛ))
    display("‖ℛ̃‖_{Z'} ≤ ")
    display(sup(ℛ̃))
    all(x -> isguaranteed(x) && isfinite(inf(x)) && isfinite(sup(x)), [wRΩ, wRΨin, RΨband, wRΨband, rRΩ̃, rRΨ̃, ℛ, ℛ̃, eN, eR, nΩ̃, rΩ̃, sΩ, sΩ̃, sΩfar, sΩ̃far, IoΩ, IoΩ̃]) || error("the bounds are not guaranteed finite intervals")
    return (wRΩ = wRΩ, wRΨin = wRΨin, RΨband = RΨband, wRΨband = wRΨband, rRΩ̃ = rRΩ̃, rRΨ̃ = rRΨ̃, nΩ̃ = nΩ̃, rΩ̃ = rΩ̃,
            sΩ = sΩ, sΩ̃ = sΩ̃, sΩfar = sΩfar, sΩ̃far = sΩ̃far, IoΩ = IoΩ, IoΩ̃ = IoΩ̃, ℛ = ℛ, ℛ̃ = ℛ̃)
end

# Check each target of Table certi1 as well as (res:final). An interval around an upper majorant is used
# only through its upper endpoint. The lower endpoint is meaningful only for the normalization.
function verify_bounds(b)
    targets = (wRΩ = I"4.3e-12", wRΨin = I"6.8e-14", RΨband = I"5.6e-13", wRΨband = I"1.25e-9",
               rRΩ̃ = I"1.25e-12", rRΨ̃ = I"1.6e-14", sΩ = I"13", sΩ̃ = I"0.025",
               sΩfar = I"5e-4", sΩ̃far = I"1e-18", IoΩ = I"0.38", IoΩ̃ = I"3.2e-5",
               ℛ = I"5.85e-12", ℛ̃ = I"1.46e-12")
    for name in keys(targets)
        x = getproperty(b, name)
        isguaranteed(x) && isfinite(sup(x)) && sup(x) ≤ inf(getproperty(targets, name)) || error("the target $name is not certified: $x")
    end
    isguaranteed(b.nΩ̃) && inf(b.nΩ̃) ≥ inf(I"1") && sup(b.nΩ̃) ≤ inf(I"1.01") || error("the normalization is not certified in [1,1.01]")
    display("All residual and coefficient targets of Table certi1 and (res:final) are certified.")
    return true
end
