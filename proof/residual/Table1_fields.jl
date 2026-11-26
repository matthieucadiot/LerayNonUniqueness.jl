# Table1_fields.jl : the raw drift and the joint test-field integral of Table 1.
# The stored odd fields are divided by z. No division by r or z is used below.

function table1_first_derivatives(Π,hu,hv,N)
    return (f=Model(Π,N),u=Model(differentiate(Π,(1,0))/hu,N),v=Model(differentiate(Π,(0,1))/hv,N))
end

function table1_glued_first(F,exterior,χ,XTL,A,B)
    E=(f=far_evaluate(exterior.f,XTL),
       u=(one(A.Π[(0,0)])/I"6")*A^3*far_evaluate(exterior.r1,XTL),
       v=(one(B.Π[(0,0)])/I"6")*B^3*far_evaluate(exterior.z1,XTL))
    return (f=E.f+χ.f*(F.f-E.f),
            u=E.u+χ.u*(F.f-E.f)+χ.f*(F.u-E.u),
            v=E.v+χ.v*(F.f-E.f)+χ.f*(F.v-E.v))
end

# Only Ψ5 and Ψ̃5 need first derivatives here; Ω* and Ω̃* need values only.
# The gluing formulas are the first three entries of glued_field, without
# constructing unused second derivatives or their far-field monomials.
function table1_profile_fields(data,far,C,R,par)
    interior=C.interior&&R.interior;N=interior ? par.N_in : par.N_band
    Ψ=table1_first_derivatives(piece(data.Ψ,C.B,C.k,R.B,R.k,data.p),C.h,R.h,N)
    Ω=(f=Model(piece(data.Ω,C.B,C.k,R.B,R.k,data.p),N),)
    Ψ̃=table1_first_derivatives(piece(data.Ψ̃,C.B,C.k,R.B,R.k,data.p),C.h,R.h,N)
    Ω̃=(f=Model(piece(data.Ω̃,C.B,C.k,R.B,R.k,data.p),N),)
    if !interior
        X,T,L=cell_variables(C.poly.A,R.poly.B,par.N_band)
        forms=(far.P.f,far.P.r1,far.P.z1,far.O.f,far.P̃.f,far.P̃.r1,far.P̃.z1,far.Õ.f)
        K=maximum(order(space(F))[1] for F in forms)
        H=maximum(order(space(F))[2] for F in forms)
        J=maximum(order(space(F))[3] for F in forms)
        XTL=monomials(X,T,L,K,H,J)
        χ=cutoff(C.χ,R.χ,C.h,R.h)
        Ψ=table1_glued_first(Ψ,far.P,χ,XTL,C.poly.A,R.poly.B)
        Ψ̃=table1_glued_first(Ψ̃,far.P̃,χ,XTL,C.poly.A,R.poly.B)
        O=far_evaluate(far.O.f,XTL);Õ=far_evaluate(far.Õ.f,XTL)
        Ω=(f=O+χ.f*(Ω.f-O),);Ω̃=(f=Õ+χ.f*(Ω̃.f-Õ),)
    end
    return (Ψ=Ψ,Ω=Ω,Ψ̃=Ψ̃,Ω̃=Ω̃,interior=interior)
end

# A³B³ times a5, b5_r/r and b5_z/z, where A=1+u and B=1+v.
# Pz/B³ = ∂z(zP), Pr/A³ = (r∂r(zP)+2zP)/z, and Qr/A³ = r∂rQ+2Q.
function table1_numerators(Ψ,Ω,Ψ̃,Ω̃,A,B)
    r2=(A^4-I"1")/I"12";z2=(B^4-I"1")/I"12"
    A3=A^3;B3=B^3
    Pz=B3*Ψ.f+I"6"*z2*Ψ.v
    Pr=I"2"*A3*Ψ.f+I"6"*r2*Ψ.u
    Qr=I"2"*A3*Ψ̃.f+I"6"*r2*Ψ̃.u
    a=A3*B3*Ω̃.f+Pz*Qr-I"6"*z2*Pr*Ψ̃.v
    br=A3*B3*Ω̃.f/I"2"-I"6"*A3*z2*Ω.f*Ψ̃.v
    bz=A3*B3*Ω̃.f/I"2"+B3*Ω.f*Qr
    return (a=a,br=br,bz=bz,Pz=Pz,Pr=Pr)
end

# [∫cell(a5²+|b5|²), sup_cell |(-r Ψ5_z,r Ψ5_r+2Ψ5)|²].
# Multiplying the three cleared numerators by A⁴B⁴ permits reuse of the residual
# densities A^-14 B^-14 J. The factors r² and z² in b5 use the corresponding densities.
function table1_cell(C,R,P,par)
    N=P.interior ? par.N_in : par.N_band
    A=response_order(C.poly.A,N);B=response_order(R.poly.B,N)
    F=table1_numerators(P.Ψ,P.Ω,P.Ψ̃,P.Ω̃,A,B)
    scale=A^4*B^4
    a=response_chart(scale*F.a,R)
    br=response_chart(scale*F.br,R)
    bz=response_chart(scale*F.bz,R)
    joint2=C.hf*R.jac*(response_integral(a,C.ρres,R.μres_even)+
                         response_integral(br,C.ρresr,R.μres_even)+
                         response_integral(bz,C.ρres,R.μres_odd))
    # The reciprocal is needed only for a pointwise model; its analytic tail is
    # enclosed by the same binomial arithmetic as all other local reciprocals.
    A=response_order(A,par.N_functional);B=response_order(B,par.N_functional)
    D=F.Pz*(B^I"-3");H=F.Pr*(A^I"-3")
    drift2=sup_bound((A^4-I"1")/I"12"*(D*D)+(B^4-I"1")/I"12"*(H*H))
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&inf(x)≥0,(joint2,drift2)) || error("invalid Table 1 cell bounds")
    return [joint2,drift2]
end

# General monomial counterpart of exterior_sup: sup |r^a z^b F| on R≥L.
function table1_exterior_sup(F,a::Integer,b::Integer,L)
    a≥0&&b≥0&&inf(L)>1 || error("invalid Table 1 exterior supremum")
    S=big(I"0")
    for (k,h,j) in indices(space(F))
        c=F[(k,h,j)]
        if !(inf(c)==0&&sup(c)==0)
            power=a+b-h
            power≤0&&sup(interval(BigFloat,power)*log(L)+interval(BigFloat,j))≤0 || error("a Table 1 exterior term is not decreasing")
            m=2*k+a
            angular=(m==0||b==0) ? big(I"1") :
                (interval(BigFloat,m)/interval(BigFloat,m+b))^(interval(BigFloat,m)/I"2")*
                (interval(BigFloat,b)/interval(BigFloat,m+b))^(interval(BigFloat,b)/I"2")
            S=S+abs(c)*angular*L^interval(BigFloat,power)*log(L)^interval(BigFloat,j)
        end
    end
    return S
end

# Exact finite exterior expressions. In particular both components of b5 and
# every term of the even potential are retained before squaring and integrating.
function table1_exterior(far,L)
    X,T,ℓ=variables()
    P=far.P.f;O=far.O.f;Q=far.P̃.f;V=far.Õ.f
    Pz=P+(I"1"-X)*𝒜z(P)
    Pr=I"2"*P+X*𝒜r(P)
    Qr=I"2"*Q+X*𝒜r(Q)
    a=V+Pz*Qr-(I"1"-X)*Pr*𝒜z(Q)
    br=V/I"2"-(I"1"-X)*O*𝒜z(Q)
    bz=V/I"2"+O*Qr
    joint2=exterior_integral(a,0,0,L)+exterior_integral(br,2,0,L)+exterior_integral(bz,0,2,L)
    drift2=table1_exterior_sup(Pz,1,0,L)^2+table1_exterior_sup(Pr,0,1,L)^2
    # The spherical exterior contains the exterior of the computational square.
    # These are upper majorants, so no positive lower endpoint is used.
    return (joint2=I"[0,1]"*joint2,drift2=drift2)
end

function assemble_table1(S,c,exterior;complete=false)
    size(S,1)==2&&size(S,2)>0 || error("Table 1 needs two quantities for every profile column")
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&inf(x)≥0,S) || error("invalid Table 1 columns")
    drift=sqrt(max(maximum(S[2,:]),exterior.drift2))
    raw=sqrt(sum(S[1,:])+exterior.joint2)
    # Ψ*=Ψ5-Δ^-1 RΨ and Ψ̃*=Ψ̃5-Δ^-1 RΨ̃, using the two estimates
    # in Section 6 for a*-a5 and b*-b5.
    aerror=c.Cb*c.RΨZ*c.nΩ̃Z+drift*c.rRΨ̃
    berror=c.sΩ*c.rRΨ̃
    correction=sqrt(aerror*aerror+berror*berror)
    test=raw+correction
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&inf(x)≥0,(drift,raw,correction,test)) || error("invalid Table 1 assembled bounds")
    return (drift=drift,raw_test=raw,correction=correction,test=test,profile_constants=c,
            certified=complete===true&&get(c,:certified,false)===true)
end

# The residual entries of Table 1 are checked separately by verify_bounds.
function verify_table1(b)
    b.certified || error("Table 1 requires complete profile columns")
    c=b.profile_constants
    for (name,x,target) in (("drift",b.drift,big(I"22.45")),("test fields",b.test,big(I"14.4")),
                          ("stream defect",c.RΨZ,big(I"1.251e-9")),("even potential gradient",c.gΨ̃,big(I"0.89")),
                          ("weighted even profile",c.RΩ̃,big(I"22.55")))
        isguaranteed(x)&&isfinite(sup(x))&&sup(x)≤inf(target) || error("Table 1 target $name is not certified: $x")
    end
    inf(c.gΨ̃_min)≥sup(big(I"0.88")) || error("Table 1 lower even-potential gradient bound is not certified")
    return true
end
