# The shared profile fields and additional Table 1 integrals.
# These routines are reused from the existing Table 2 code.

function profile_fields_on_cell(data,far,C,R,par)
    interior=C.interior&&R.interior;N=interior ? par.N_in : par.N_band
    Ψ=field(piece(data.Ψ,C.B,C.k,R.B,R.k,data.p),C.h,R.h,N)
    Ω=field(piece(data.Ω,C.B,C.k,R.B,R.k,data.p),C.h,R.h,N)
    Ψ̃=field(piece(data.Ψ̃,C.B,C.k,R.B,R.k,data.p),C.h,R.h,N)
    Ω̃=field(piece(data.Ω̃,C.B,C.k,R.B,R.k,data.p),C.h,R.h,N)
    χ=cutoff(C.χ,R.χ,C.h,R.h);XTL=nothing
    if !interior
        X,T,L=cell_variables(C.poly.A,R.poly.B,par.N_band)
        XTL=monomials(X,T,L,far.K,far.H,far.J)
        Ψ=glued_field(Ψ,far_field_on_cell(far.P,XTL,C.poly.A,R.poly.B),χ)
        Ω=glued_field(Ω,far_field_on_cell(far.O,XTL,C.poly.A,R.poly.B),χ)
        Ψ̃=glued_field(Ψ̃,far_field_on_cell(far.P̃,XTL,C.poly.A,R.poly.B),χ)
        Ω̃=glued_field(Ω̃,far_field_on_cell(far.Õ,XTL,C.poly.A,R.poly.B),χ)
    end
    return (Ψ=Ψ,Ω=Ω,Ψ̃=Ψ̃,Ω̃=Ω̃,χ=χ,interior=interior,XTL=XTL,λ̃=data.λ̃)
end

response_chart(M,R)=R.axis ? substitute_axis(M) : M
response_order(M,N)=Model(M.Π,M.ϱ,N)

function response_gradient(G,A,B,odd::Bool)
    r1=I"6"*(A^I"-3")*G.u
    z1=I"6"*(B^I"-3")*G.v
    if odd
        z1=G.f+((B^4-I"1")/I"12")*z1
    end
    return (r1=r1,z1=z1)
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

# All quantities on a cell share the same four profile fields.
function profile_cell(C,R,data,far,par)
    P=profile_fields_on_cell(data,far,C,R,par)
    return vcat(cell_bounds(C,R,P,data,par),profile_extra_cell(C,R,P,par),table1_cell(C,R,P,par))
end
