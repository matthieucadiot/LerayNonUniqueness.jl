# The additional scalar integrals; all use the same density/moment machinery as the profile residuals.

function response_farfield(j,data,far)
    data.odd && return nothing
    g=response_multipoles(data.a[j,:])
    return (g=g,jet=far_field_derivatives(g),ext=response_exterior(g,far.O.f,far.Õ.f,data.μ[j]))
end

function profile_extra_cell(C,R,P,par)
    A=response_order(C.poly.A,par.N_ρ);B=response_order(R.poly.B,par.N_ρ)
    N=response_chart(numerator_Ψ(P.Ψ,P.Ω,C.poly,R.poly),R)
    grad=response_gradient(P.Ψ̃,A,B,false)
    gr=response_chart(grad.r1,R);gz=response_chart(grad.z1,R);Ω=response_chart(P.Ω̃.f,R)
    jac=C.hf*R.jac
    nRΨ=jac*response_integral(N,C.ρres,R.μres_odd)
    gΨ̃=jac*(response_integral(gr,C.ρJr,R.μJeven)+response_integral(gz,C.ρJ,R.μJodd))
    RΩ̃=jac*(response_integral(Ω,C.ρJr,R.μJeven)+response_integral(Ω,C.ρJ,R.μJodd))
    return [nRΨ,gΨ̃,RΩ̃]
end

function profile_extra_exterior(far,L)
    g=far.P̃
    return (gΨ̃2_ext=exterior_integral(g.r1,2,0,L)+exterior_integral(g.z1,0,2,L),
            RΩ̃2_ext=exterior_integral(far.Õ.f,2,0,L)+exterior_integral(far.Õ.f,0,2,L))
end

# At u=0 the square root is not analytic in the affine u coordinate. A range model is sufficient there;
# functional_models handles the support boundary without dividing by the radial coordinate at zero.
function functional_coordinates(C,R,par)
    N=par.N_functional
    A=response_order(C.poly.A,N)
    U=(A^4-I"1")/I"12"
    Ir=interval(inf(sqrt(max(I"0",U0(C.a)))),sup(sqrt(U0(C.b))))
    Iz=interval(inf(sqrt(max(I"0",U0(R.c)))),sup(sqrt(U0(R.e))))
    r=inf(C.U0a)>0 ? U^I"0.5" : functional_range_model(Ir,A)
    z=functional_order(R.z,N)
    return r,z,Ir,Iz
end

function functional_cell_integrals(Q,odd,C,R)
    jac=C.hf*R.jac
    q2=jac*response_integral(Q.q,C.ρJ,R.μJeven)
    if odd
        z=R.zsup;r2=C.U0b
        k2=z*(I"1"+z)*(I"1"+z)*(I"2"+z)/I"6"
        hardy=r2*k2/(r2+k2)*q2
    else
        hardy=jac*response_integral(Q.q,C.ρJr,R.μJeven)
    end
    vector=jac*(response_integral(Q.Vr,C.ρJ,R.μJeven)+response_integral(Q.Vz,C.ρJ,R.μJeven))
    return [hardy,vector]
end

# Integrate the scalar primitive only on cells that are not wholly beyond its radial support.
# The remaining entire half-cylinder r≥upper(b) is added analytically by functional_tail_integral.
function functional_on_cell(f,bases,C,R,par)
    Ir=interval(inf(sqrt(max(I"0",U0(C.a)))),sup(sqrt(U0(C.b))))
    Iz=interval(inf(sqrt(max(I"0",U0(R.c)))),sup(sqrt(U0(R.e))))
    outside= sup(Ir)≤inf(f.r_spec[1]) || inf(Ir)≥sup(f.r_spec[4]) ||
             sup(Iz)≤inf(f.z_spec[1]) || inf(Iz)≥sup(f.z_spec[4])
    if outside
        z=functional_constant(I"0",response_order(C.poly.A,par.N_functional))
        return (q=z,qnum=z,δq=z,Vr=z,Vz=z,tail=z)
    end
    r,z,Ir,Iz=functional_coordinates(C,R,par)
    return functional_models(f,bases,r,z,Ir,Iz;N=par.N_functional)
end

# The only noncompact functional primitive is scalar radial mode m=0: Vr=C(z)/r³.
# L is chosen as the upper endpoint of the radial support, not the physical square's radius.
# Its entire remaining radial half-cylinder is integrated exactly, with radial integral 1/(2L²).
function functional_tail_integral(f,bases,row,par,L)
    L=interval(Float64,L)
    f.kind=="scalar" || return I"0"
    inf(L)≥sup(f.r_spec[4]) || error("the scalar primitive tail radius must be beyond its support")
    !isempty(row) && sup(row[1].c)==0 && inf(row[end].zsup)≥sup(f.z_spec[4]) || error("the scalar tail rows do not cover the positive vertical support")
    all(isequal_interval(row[j].e,row[j+1].c) for j in 1:length(row)-1) || error("the scalar tail row partition has a gap")
    mass=I"0"
    for R in row
        N=par.N_functional
        r=functional_constant(I"2"*L,response_order(R.poly.B,N));z=functional_order(R.z,N)
        Ir=I"2"*L;Iz=interval(inf(sqrt(max(I"0",R.U0c))),sup(R.zsup))
        Q=functional_models(f,bases,r,z,Ir,Iz;N=N)
        # x integration contributes2; dz=(1/6) B³/sqrt(U0(v)) dv.
        unit=Model(Sequence(Taylor(0) ⊗ Taylor(0),[I"1"]),N)
        ρ=density(unit,par.dmom)
        mass=mass+R.jac*response_integral(Q.tail,ρ,R.μJeven)/I"12"
    end
    return I"2"*(I"1"*π)^2/(L*L)*mass
end

function raw_response_bounds(S,odd,ext;complete=false)
    if odd
        return (wE=sqrt(S[1]),wD_in=sqrt(S[2]),D_band=sqrt(S[3]),wD_band=sqrt(S[4]),nΦ=sqrt(S[5]),gΨ=sqrt(S[6]),complete=complete)
    end
    return (rE=sqrt(S[1]+I"[0,1]"*ext.rE2),rD=sqrt(S[2]),nΦ=sqrt(S[3]),RΦ=sqrt(S[4]),
            gΨ=sqrt(S[5]+I"[0,1]"*ext.gPsi2),orth=S[6],complete=complete)
end
