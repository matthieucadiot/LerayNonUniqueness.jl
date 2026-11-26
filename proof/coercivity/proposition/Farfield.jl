# Proposition 4.21 (i): the coefficients of the raw potential Ψ₅, on the
# finite computational square and on its infinite exterior. The square still
# contains spline and gluing cells with R≥200; the exterior bound alone does
# not cover those cells.

# If Ψ=zP(x,t,ℓ), then ΔΨ=z𝓛(P,1). The fourth coefficient therefore comes
# from ΔΨ, rather than from the independently stored approximation Ω₅.
function coercivity_series(P)
    x,t,ℓ=variables();y=I"1"-x
    Ar=𝒜r(P);Az=𝒜z(P)
    rz=x*(Ar+y*(𝒜r(Az)-I"2"*Az))
    z=P+y*Az
    O=𝓛(P,1)
    return (err=-rz-I"2"*z,ezz=rz+z,
            A12=(I"4"*Ar+x*(𝒜r(Ar)-I"2"*Ar)-I"3"*Az-y*(𝒜z(Az)-I"2"*Az))/I"2",
            e0=O+y*𝒜z(O))
end

# sup_{R≥R₀} R^shift |F|, with x∈[0,1]. Each radial monomial is
# decreasing on [R₀,∞); the logarithmic factors are retained explicitly.
function farfield_series_bound(F,R₀,shift)
    inf(R₀)>1 || error("the exterior radius must exceed 1")
    all(isguaranteed,coefficients(F)) || error("a far-field coefficient lost its interval guarantee")
    B=I"0"
    for (k,h,j) in indices(space(F))
        c=F[(k,h,j)]
        inf(c)==0 && sup(c)==0 && continue
        h>shift || error("the scaled far-field tail does not decay")
        sup(interval(Float64,j))≤inf(interval(Float64,h-shift)*log(R₀)) || error("a logarithmic tail is not decreasing at the exterior radius")
        B=B+abs(c)*R₀^interval(Float64,shift-h)*log(R₀)^interval(Float64,j)
    end
    return B
end

# This transcription check uses Ψ₋₁=130z/R² and the four exact formulas
# (4.65). A₁₂ is the coefficient after removing sqrt(x(1-x)).
function check_farfield_series()
    x,t,ℓ=variables();y=I"1"-x;σ=I"130"
    F=coercivity_series(σ*t^2*(I"1"+I"0"*x))
    lead=(err=I"2"*σ*t^2*y*(I"4"*y-I"3"),
          ezz=σ*t^2*(-x*x+I"6"*x*y-y*y),
          A12=σ*t^2*(I"8"*x-I"5"),
          e0=-I"6"*σ*t^4*(x-I"3"*y))
    for name in keys(lead)
        all(c->in_interval(0,c),coefficients(getproperty(F,name)-getproperty(lead,name))) || error("the leading far-field formulas disagree")
    end
    return true
end

# The exact finite expansion is used only outside the computational square.
# Its coefficients and differential operators are reused from Asymptotics.jl.
# The factor sqrt(x(1-x)) of a₁₂ is bounded by 1/2.
function exterior_farfield(R₀)
    R₀=interval(Float64,R₀);σ=I"130"
    check_farfield_series()
    x,t,ℓ=variables()
    P=far_field(precision(Float64)).P
    F=coercivity_series(P-σ*t^2)
    B=(err=farfield_series_bound(F.err,R₀,2)/σ,
       ezz=farfield_series_bound(F.ezz,R₀,2)/σ,
       a12=farfield_series_bound(F.A12,R₀,2)/(I"2"*σ),
       e0=farfield_series_bound(F.e0,R₀,4)/σ)
    ok=sup(B.err)≤inf(I"0.05") && sup(B.ezz)≤inf(I"0.05") &&
       sup(B.a12)≤inf(I"0.05") && sup(B.e0)≤1
    return (ok=ok,bounds=B,radius=R₀)
end

# On a spline or gluing cell, compare the full coefficients with their
# leading terms before taking absolute values. Bounds are normalized by σ.
function farfield_box_bounds(F)
    σ=I"130";R²=functional_range(F.R2)
    invR⁶=F.R2^I"-3"
    lead=(err=I"2"*σ*F.z2*(F.z2-I"3"*F.r2)*invR⁶,
          ezz=σ*(-F.r2*F.r2+I"6"*F.r2*F.z2-F.z2*F.z2)*invR⁶,
          a12=-σ*F.r*F.z*(I"5"*F.z2-I"3"*F.r2)*invR⁶,
          e0=-I"6"*σ*(F.r2-I"3"*F.z2)*invR⁶)
    return (err=abs(functional_range(F.err-lead.err))*R²/σ,
            ezz=abs(functional_range(F.ezz-lead.ezz))*R²/σ,
            a12=abs(functional_range(F.a12-lead.a12))*R²/σ,
            e0=abs(functional_range(F.e0-lead.e0))*R²*R²/σ)
end

function check_farfield_box(F)
    sup(functional_range(F.R2))<40000 && return true
    B=farfield_box_bounds(F)
    return sup(B.err)≤inf(I"0.05") && sup(B.ezz)≤inf(I"0.05") &&
           sup(B.a12)≤inf(I"0.05") && sup(B.e0)≤1
end
