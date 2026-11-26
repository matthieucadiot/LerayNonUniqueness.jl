# Response_norms.jl : scalar cell quantities for (6.19), (6.21), (6.22), using the residual integration code.

# Cache the four profile fields once for all responses on a cell. Derivatives are always in u,v.
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

# Signed pairing: no truncation of the polynomial product and no positivity intersection.
# |f g−Π_f Π_g|≤ε_Q; |ρ μ−Π_ρ Π_μ|≤δ; ∫|Π_f Π_g|≤4|Π_f Π_g|_1.
function cell_pairing(F::Model,G::Model,ρ,μ)
    Π=F.Π*G.Π
    n1,n2=order(space(Π))
    length(ρ.m)>n1 && length(μ.m)>n2 || error("not enough moments for the signed pairing")
    value=I"0"
    for (i,j) in indices(space(Π))
        value=value+Π[(i,j)]*ρ.m[i+1]*μ.m[j+1]
    end
    ε=F.ϱ*norm(G.Π,1)+G.ϱ*norm(F.Π,1)+F.ϱ*G.ϱ
    δ=ρ.ϱ*μ.Π1+μ.ϱ*ρ.Π1+ρ.ϱ*μ.ϱ
    err=I"4"*δ*norm(Π,1)+ε*(ρ.m[1]*μ.m[1]+I"4"*δ)
    return value+I"[-1,1]"*err
end

response_chart(M,R)=R.axis ? substitute_axis(M) : M
response_order(M,N)=Model(M.Π,M.ϱ,N)

# The cutoff is applied after the rational weight. Even potentials are glued to the harmonic multipoles.
function response_cell_fields(j,rdata,CR,RR,C,R,P,farresponse,par)
    N=P.interior ? par.N_response : par.N_band
    A=response_order(C.poly.A,N);B=response_order(R.poly.B,N)
    P_G,P_F=response_spline_jets(j,rdata,CR,RR,par.N_in)
    G=weighted_response_jets(P_G,A,B,rdata.alpha_G)
    F=weighted_response_jets(P_F,A,B,rdata.alpha_F)
    if !P.interior
        F=response_product(F,P.χ)
        if rdata.odd
            G=response_product(G,P.χ)
        else
            farresponse===nothing && error("an even band potential needs its exterior multipoles")
            ext=far_field_on_cell(farresponse.jet,P.XTL,C.poly.A,R.poly.B)
            G=glued_field(G,ext,P.χ)
        end
    end
    return (G=G,F=F,P_G=P_G,P_F=P_F,A=A,B=B,N=N)
end

# Odd output: [||sqrt(w) E||², ||sqrt(w) Eψ_in||², ||Eψ_band||²,
#              ||sqrt(w) Eψ_band||², ||Φ||², ||∇Ψ||²] on this cell.
# Even output: [||r E||², ||r Eψ||², ||Φ||², ||R Φ||², ||∇Ψ||², <Φ,Ψ̃_5>].
# q is the PHYSICAL source model, already in the axis chart when R.axis is true.
function response_cell_quantities(j,rdata,CR,RR,C,R,P,farresponse,q,par)
    precision(BigFloat)==par.input_prec || error("response cell precision differs from the working precision")
    W=response_cell_fields(j,rdata,CR,RR,C,R,P,farresponse,par)
    F=W.F;G=W.G;a=C.poly;b=R.poly
    zero=I"0"*P.Ω.f
    if P.interior
        if rdata.odd
            main=cleared_response_numerator_odd(P.Ψ,P.Ω,W.P_G,W.P_F,zero,a,b)
        else
            main=cleared_response_numerator_even(P.Ψ,P.Ω,P.Ω̃,W.P_G,W.P_F,zero,P.λ̃,rdata.μ[j],a,b)
        end
        defect=cleared_response_defect(W.P_G,W.P_F,a,b,rdata.odd)
        E=main.N*(response_order(main.H,W.N)^(-interval(Float64,main.power)))
        Eψ=defect.N*(response_order(defect.H,W.N)^(-interval(Float64,defect.power)))
    elseif rdata.odd
        E=response_numerator_odd(P.Ψ,P.Ω,G,F,zero,a,b)
        Eψ=numerator_Ψ(G,F,a,b)
    else
        E=response_numerator_even(P.Ψ,P.Ω,P.Ω̃,G,F,zero,P.λ̃,rdata.μ[j],a,b)
        Eψ=numerator_Ψ̃(G,F,a,b)
    end
    gradient=response_gradient(G,W.A,W.B,rdata.odd)
    E=response_chart(E,R);Eψ=response_chart(Eψ,R)
    D=response_chart(a.D*b.D,R)
    if rdata.odd
        E=R.z*E                                                     # the odd equation's physical z factor is inserted before subtracting q
    end
    E=E-D*q
    Φ=response_chart(F.f,R)
    gr=response_chart(gradient.r1,R);gz=response_chart(gradient.z1,R)
    jac=C.hf*R.jac
    if rdata.odd
        z=R.zsup;r2=C.U0b
        k2=z*(I"1"+z)*(I"1"+z)*(I"2"+z)/I"6"
        w=r2*k2/(r2+k2)
        e=jac*w*response_integral(E,C.ρres,R.μres_even)
        d=jac*response_integral(Eψ,C.ρres,R.μres_odd)
        n=jac*response_integral(Φ,C.ρJ,R.μJodd)
        g=jac*(response_integral(gr,C.ρJr,R.μJodd)+response_integral(gz,C.ρJ,R.μJeven))
        o=I"0"
        values=[e,P.interior ? w*d : o,P.interior ? o : d,P.interior ? o : w*d,n,g]
    else
        e=jac*response_integral(E,C.ρresr,R.μres_even)
        d=jac*response_integral(Eψ,C.ρresr,R.μres_even)
        n=jac*response_integral(Φ,C.ρJ,R.μJeven)
        r=jac*(response_integral(Φ,C.ρJr,R.μJeven)+response_integral(Φ,C.ρJ,R.μJodd))
        g=jac*(response_integral(gr,C.ρJr,R.μJeven)+response_integral(gz,C.ρJ,R.μJodd))
        pair=jac*response_pairing(Φ,response_chart(P.Ψ̃.f,R),C.ρJ,R.μJeven)
        values=[e,d,n,r,g,pair]
    end
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),values) || error("a response cell quantity is not guaranteed finite")
    all(x->precision(inf(x))==precision(sup(x))==par.prec,values) || error("a response cell quantity has the wrong precision")
    return values
end

function exterior_response_quantities(farresponse,L)
    ext=farresponse.ext
    return (rE2=exterior_integral(ext.E,2,0,L),
            gPsi2=exterior_integral(ext.r1,2,0,L)+exterior_integral(ext.z1,0,2,L))
end
