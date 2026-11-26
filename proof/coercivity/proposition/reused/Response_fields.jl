# Response_fields.jl : response fields and local errors of (6.18), with the representation (6.24).
# Odd fields are divided by z throughout; even fields have no parity factor. All derivative entries are in u,v.

function response_product(F,G)
    return (f=F.f*G.f,
            u=F.u*G.f+F.f*G.u, v=F.v*G.f+F.f*G.v,
            uu=F.uu*G.f+I"2"*F.u*G.u+F.f*G.uu,
            vv=F.vv*G.f+I"2"*F.v*G.v+F.f*G.vv)
end

function response_zero(F)
    z=I"0"*F.f
    return (f=z,u=z,v=z,uu=z,vv=z)
end

# H=1+R², H_u=A³/3, H_uu=A², with the corresponding v derivatives.
function response_geometry(A,B)
    H=I"1"+(A^4+B^4-I"2")/I"12"
    return (H=H,Hu=A^3/I"3",Hv=B^3/I"3",Huu=A^2,Hvv=B^2)
end

# The analytic weight H^-s and its derivatives are enclosed separately. A bounded model remainder is never differentiated.
function response_weight(A,B,s::Integer)
    s≥0 || error("the response weight exponent must be nonnegative")
    g=response_geometry(A,B); H=g.H
    α=-interval(Float64,s)
    W0=H^α; W1=H^(α-I"1"); W2=H^(α-I"2")
    return (f=W0,u=α*W1*g.Hu,v=α*W1*g.Hv,
            uu=α*(α-I"1")*W2*g.Hu*g.Hu+α*W1*g.Huu,
            vv=α*(α-I"1")*W2*g.Hv*g.Hv+α*W1*g.Hvv)
end

weighted_response_field(Π,hu,hv,A,B,s,N)=response_product(field(Π,hu,hv,N),response_weight(A,B,s))
weighted_response_jets(P,A,B,s)=response_product(P,response_weight(A,B,s))

function response_power(H,n::Integer)
    n≥0 || error("a cleared response power must be nonnegative")
    return n==0 ? Model(one(H.Π),H.N) : H^n
end
clear_response_jets(P,A,B,s,d)=cleared_response_field(P,A,B,s,d)

# H^d times the five derivatives of P/H^s. For d≥s+2 all five entries are polynomial models.
# The entries are scaled derivatives, not the derivatives of the scaled f entry.
function cleared_response_field(P,A,B,s::Integer,d::Integer)
    s≥0 && d≥s+2 || error("cleared response jets require s≥0 and d≥s+2")
    g=response_geometry(A,B); H=g.H; c=interval(Float64,s)
    return (f=response_power(H,d-s)*P.f,
            u=response_power(H,d-s-1)*(H*P.u-c*g.Hu*P.f),
            v=response_power(H,d-s-1)*(H*P.v-c*g.Hv*P.f),
            uu=response_power(H,d-s-2)*(H*H*P.uu-I"2"*c*H*g.Hu*P.u+(c*(c+I"1")*g.Hu*g.Hu-c*H*g.Huu)*P.f),
            vv=response_power(H,d-s-2)*(H*H*P.vv-I"2"*c*H*g.Hv*P.v+(c*(c+I"1")*g.Hv*g.Hv-c*H*g.Hvv)*P.f))
end

# D N5(zP,zO)/z = 2D N[P,O] of (u:N). It is linear in each input.
function response_coupling_odd(P,O,a,b)
    return I"4"*a.K*b.K*(P.u*O.v-P.v*O.u)-I"2"*a.K*b.D*(P.f*O.u-P.u*O.f)+
           I"4"*a.D*b.K*P.f*O.v+I"2"*a.D*b.D*P.f*O.f
end

# D E_j/z and D Ẽ_j of (6.18). q is the computational functional (q/z for odd parity).
function response_numerator_odd(Ψ,Ω,G,F,q,a,b)
    return numerator_Ω(Ψ,F,a,b)+response_coupling_odd(G,Ω,a,b)-a.D*b.D*q
end
function response_numerator_even(Ψ,Ω,Ω̃,G,F,q,λ̃,μ̃,a,b)
    return numerator_Ω̃(Ψ,Ω,G,F,λ̃,a,b)-a.D*b.D*(μ̃*Ω̃.f+q)
end

# Preserve cancellation in the interior: clear H before adding all terms, then apply the common reciprocal.
# P_F,P_G are the unweighted polynomial jets of the STORED weighted spline numerators.
function cleared_response_numerator_odd(Ψ,Ω,P_G,P_F,q,a,b)
    H=response_geometry(a.A,b.B).H; d=6
    F=cleared_response_field(P_F,a.A,b.B,4,d)
    G=cleared_response_field(P_G,a.A,b.B,2,d)
    return (N=response_numerator_odd(Ψ,Ω,G,F,response_power(H,d)*q,a,b),H=H,power=d)
end
function cleared_response_numerator_even(Ψ,Ω,Ω̃,P_G,P_F,q,λ̃,μ̃,a,b)
    H=response_geometry(a.A,b.B).H; d=5
    F=cleared_response_field(P_F,a.A,b.B,3,d)
    G=cleared_response_field(P_G,a.A,b.B,1,d)
    N=numerator_Ω̃(Ψ,Ω,G,F,λ̃,a,b)-a.D*b.D*response_power(H,d)*(μ̃*Ω̃.f+q)
    return (N=N,H=H,power=d)
end

function cleared_response_defect(P_G,P_F,a,b,odd::Bool)
    sG=odd ? 2 : 1; sF=odd ? 4 : 3; d=odd ? 4 : 3
    H=response_geometry(a.A,b.B).H
    G=cleared_response_field(P_G,a.A,b.B,sG,d)
    F=(f=response_power(H,d-sF)*P_F.f,)
    N=odd ? numerator_Ψ(G,F,a,b) : numerator_Ψ̃(G,F,a,b)
    return (N=N,H=H,power=d)
end

# Physical gradients of the potential. For odd zG the radial term is r*z*r1 and axial term z1.
# For even G the radial term is r*r1 and the axial term is z*z1.
function response_gradient(G,A,B,odd::Bool)
    r1=I"6"*(A^I"-3")*G.u
    z1=I"6"*(B^I"-3")*G.v
    if odd
        z1=G.f+((B^4-I"1")/I"12")*z1
    end
    return (r1=r1,z1=z1)
end

# Harmonic multipoles (6.23): a_m h_2m R^(-4m-3)=a_m t^(2m+3) p_m(x).
function response_multipoles(a)
    length(a)==4 && all(isguaranteed,a) || error("four guaranteed interval multipole coefficients are required")
    a=interval.(Float64,a)
    X,T,ℓ=variables()
    p=(one(X),I"4"-I"5"*X,I"8"-I"28"*X+I"21"*X^2,
       I"64"-I"432"*X+I"792"*X^2-I"429"*X^3)
    return a[1]*T^3*p[1]+a[2]*T^5*p[2]+a[3]*T^7*p[3]+a[4]*T^9*p[4]
end

# The only nonzero exterior error is (6.25); the potential gradient is supplied by far_field_derivatives(g).
function response_exterior(g,O,Ω̃,μ̃)
    μ̃=interval(Float64,μ̃)
    X,T,ℓ=variables()
    Arg=𝒜r(g);Azg=𝒜z(g)
    E=(X*Arg+I"2"*g)*(O+(I"1"-X)*𝒜z(O))-X*(I"1"-X)*Azg*𝒜r(O)-μ̃*Ω̃
    return (E=E,r1=T*T*Arg,z1=T*T*Azg)
end
