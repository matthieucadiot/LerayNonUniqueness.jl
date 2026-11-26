# Coefficients (4.11): Ψ₅=zP. In particular e₀=∂zΔΨ₅; Ω₅ is not substituted for ΔΨ₅.
# Differentiate polynomial pieces and explicit analytic formulas BEFORE adding model errors.

coordinate(r) = (I"1"+I"12"*r*r)^I"0.25"-I"1"

function affine_model(a,b,axis,N)
    P=zeros(Interval{Float64},axis==1 ? Taylor(1)⊗Taylor(0) : Taylor(0)⊗Taylor(1))
    P[(0,0)]=proof_float64((a+b)/I"2")
    P[axis==1 ? (1,0) : (0,1)]=proof_float64((b-a)/I"2")
    return Model(P,N)
end

# All derivatives ∂u^i∂v^j P, i+j≤3, from the SAME local polynomial.
function polynomial_jets(P,hu,hv,N)
    return Dict((i,j)=>Model((one(first(coefficients(P)))/(hu^interval(i)*hv^interval(j)))*differentiate(P,(i,j)),N)
                for i in 0:3 for j in 0:3-i)
end

function product_jets(P,Q)
    return Dict((i,j)=>sum(interval(Float64,binomial(i,a)*binomial(j,b))*P[(a,b)]*Q[(i-a,j-b)]
                         for a in 0:i for b in 0:j) for i in 0:3 for j in 0:3-i)
end

function series_model(F,X,T,L,N,XTL=nothing)
    # The residual code already supplies the monomial arithmetic and its errors.
    # Reuse each monomial across all derivative series on the same cell.
    if XTL===nothing
        K,H,J=order(space(F))
        XTL=monomials(X,T,L,K,H,J)
    end
    return far_evaluate(F,XTL)
end

# Conversion from D_r=(1/r)∂r and D_z=(1/z)∂z to u,v derivatives, through order three.
function derivative_conversion(A)
    Z=I"0"*A; O=functional_constant(I"1",A)
    return ((O,), (Z,A^3/I"6"), (Z,A^2/I"2",A^6/I"36"),
            (Z,A,A^5/I"4",A^9/I"216"))
end

function exterior_jets(A,B,far,N)
    r2=(A^4-I"1")/I"12"; z2=(B^4-I"1")/I"12"; R2=r2+z2
    X=r2*R2^I"-1"; T=R2^I"-0.5"; L=log(R2)/I"2"
    _,t,_=variables()
    F=Dict((0,0)=>far.P)
    for i in 0:3, j in 0:3-i
        i+j==0 && continue
        F[(i,j)]= i>0 ? t*t*𝒜r(F[(i-1,j)]) : t*t*𝒜z(F[(0,j-1)])
    end
    K=maximum(order(space(v))[1] for v in values(F))
    H=maximum(order(space(v))[2] for v in values(F))
    J=maximum(order(space(v))[3] for v in values(F))
    XTL=Array{Model}(undef,K+1,H+1,J+1)
    cache=Dict((0,0,0)=>functional_constant(I"1",X))
    function monomial(k,h,j)
        return get!(cache,(k,h,j)) do
            j>0 ? monomial(k,h,j-1)*L : k>0 ? monomial(k-1,h,0)*X : monomial(0,h-1,0)*T
        end
    end
    # Build only monomials used by a nonzero coefficient in one of the ten jets.
    for G in values(F), (k,h,j) in indices(space(G))
        isthinzero(G[(k,h,j)]) || monomial(k,h,j)
    end
    for ((k,h,j),value) in cache
        XTL[k+1,h+1,j+1]=value
    end
    E=Dict(k=>series_model(v,X,T,L,N,XTL) for (k,v) in F)
    U=derivative_conversion(A); V=derivative_conversion(B)
    return Dict((i,j)=>sum(U[i+1][a+1]*V[j+1][b+1]*E[(a,b)]
                         for a in 0:i for b in 0:j
                         if !isthinzero(sup_bound(U[i+1][a+1])) && !isthinzero(sup_bound(V[j+1][b+1])))
                for i in 0:3 for j in 0:3-i)
end

function profile_jets(data,box,N,far)
    a,b,c,d=box; hu=proof_float64((b-a)/I"2"); hv=proof_float64((d-c)/I"2")
    ku=knot_interval(data.tu,a); kv=knot_interval(data.tv,c)
    Bu=local_basis(data.tu,data.p,ku,a,b); Bv=local_basis(data.tv,data.p,kv,c,d)
    P=polynomial_jets(piece(data.Ψ,Bu,ku,Bv,kv,data.p),hu,hv,N)
    A=I"1"+affine_model(a,b,1,N); B=I"1"+affine_model(c,d,2,N)
    if sup(b)>data.L-data.d || sup(d)>data.L-data.d
        χ=cutoff_u(a,b,data.L,data.d,100)*cutoff_v(c,d,data.L,data.d,100)
        isthinzero(χ.ϱ) || error("cutoff polynomial was truncated")
        C=polynomial_jets(χ.Π,hu,hv,N)
        E=exterior_jets(A,B,far,N)
        D=Dict(k=>P[k]-E[k] for k in keys(P))
        D=product_jets(C,D)
        P=Dict(k=>E[k]+D[k] for k in keys(P))
    end
    # v=s² at z=0 removes the square root singularity. Derivatives remain u,v derivatives.
    if isthinzero(c)
        P=Dict(k=>substitute_axis(v) for (k,v) in P)
        A=substitute_axis(A); B=substitute_axis(B)
        s=affine_model(I"0",sqrt(d),2,N)
        H=(I"4"+I"6"*s^2+I"4"*s^4+s^6)/I"12"
        z=s*H^I"0.5"
        dz=B^3*H^I"-0.5"/I"3"
        jac=proof_float64((b-a)*sqrt(d)/I"4")
    else
        z=((B^4-I"1")/I"12")^I"0.5"
        dz=B^3*z^I"-1"/I"6"
        jac=proof_float64((b-a)*(d-c)/I"4")
    end
    r2=(A^4-I"1")/I"12"
    Ir=proof_float64(sqrt(intersect_interval(U0(interval(inf(a),sup(b))),I"[0,inf]")))
    R=functional_range(r2)
    r=inf(R)>0 && sup(R)<4*inf(R) ? r2^I"0.5" : functional_range_model(Ir,A)
    return (P=P,A=A,B=B,r=r,z=z,r2=r2,z2=(B^4-I"1")/I"12",dz=dz,jac=jac,Ir=Ir)
end

function coefficients_from_jets(P,A,B,r,z)
    r2=(A^4-I"1")/I"12"; z2=(B^4-I"1")/I"12"
    Ar=A^I"-3"; Bz=B^I"-3"
    Ψz=P[(0,0)]+I"6"*z2*Bz*P[(0,1)]
    rΨrz=I"6"*r2*Ar*(P[(1,0)]+I"6"*z2*Bz*P[(1,1)])
    Ψrr=z*(I"36"*r2*A^I"-6"*P[(2,0)]+I"3"*(I"3"-A^4)*A^I"-7"*P[(1,0)])
    Ψzz=z*(I"36"*z2*B^I"-6"*P[(0,2)]+I"9"*(B^4+I"1")*B^I"-7"*P[(0,1)])
    Ψr=I"6"*r*z*Ar*P[(1,0)]
    cr2=I"3"*(A^4-I"1")*A^I"-6"
    cr1=I"3"*(I"5"*A^4+I"3")*A^I"-7"
    cz2=I"3"*(B^4-I"1")*B^I"-6"
    cz1=I"9"*(B^4+I"1")*B^I"-7"
    ΔP=cr2*P[(2,0)]+cr1*P[(1,0)]+cz2*P[(0,2)]+cz1*P[(0,1)]
    ΔPv=cr2*P[(2,1)]+cr1*P[(1,1)]+cz2*P[(0,3)]+
         (-I"6"*B^I"-3"+I"18"*B^I"-7"+cz1)*P[(0,2)]+
         (-I"27"*B^I"-4"-I"63"*B^I"-8")*P[(0,1)]
    return (err=-rΨrz-I"2"*Ψz,ezz=rΨrz+Ψz,a12=(r*(Ψrr-Ψzz)+I"3"*Ψr)/I"2",
            e0=ΔP+I"6"*z2*Bz*ΔPv)
end

function coefficient_models(data,box,N,far)
    C=profile_jets(data,box,N,far)
    E=coefficients_from_jets(C.P,C.A,C.B,C.r,C.z)
    return merge(C,E,(R2=C.r2+C.z2,u=C.A-I"1",v=C.B-I"1",measure=C.jac*C.r2*C.A^3*C.dz/I"6"))
end

# M_i−η_i I. Near r=0 the singular positive diagonal is replaced by a LOWER bound.
# This substitute is used only for positivity; Gram integration requires the full model.
function matrix_field(F,row,ζ=I"0"; positivity=false)
    R=functional_range(F.r2)
    invr2=if inf(R)>0 && sup(R)<4*inf(R)
        F.r2^I"-1"
    else
        positivity || error("subdivide the radial cell for the Gram integral")
        functional_constant(one(F.Ir)/interval(Float64,sup(F.Ir))^2,F.A)
    end
    if row.s==1
        return (F.e0+I"0.035"*invr2-row.η,)
    end
    c=I"0.25"-I"0.100001"
    row.block==:even && (c=c+I"0.235059685902438"-I"0.224")
    a=c+I"3"*(I"1"-I"0.100001")*invr2+F.err-row.η
    b=F.a12
    d=c+F.ezz-row.η-(row.block==:odd ? ζ : I"0")
    return (a,b,d)
end
