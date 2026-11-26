# Interior odd vector fields. Truncate only with the validated uniform Taylor
# remainder. The high-precision B-spline basis tables are enclosed in binary64
# before the two rigorous interval matrix multiplications.
const FIELD_ORDER=6
reduce_order(x)=functional_order(x,FIELD_ORDER)
function radial_constants(data,row,a,b,V)
    F=coefficient_models(data.profile,(a,b,V[1],V[2]),12,data.far)
    A,r2,r=F.A,F.r2,F.r
    @assert (sup(b)<data.profile.L-data.profile.d)
    ku=knot_interval(data.profile.tu,a)
    Bu=proof_float64(local_basis(data.profile.tu,data.profile.p,ku,a,b))
    return (;F,ku,Bu,hu=proof_float64((b-a)/I"2"),
            Rzr=reduce_order(I"6"*r2*A^I"-3"),
            Rrr1=reduce_order(I"36"*r2*A^I"-6"),
            Rrr2=reduce_order(I"3"*(I"3"-A^4)*A^I"-7"),
            Rr1=reduce_order(I"6"*r*A^I"-3"),r=reduce_order(r),
            basea=reduce_order(I"0.25"-I"0.100001"+I"3"*(I"1"-I"0.100001")*(r2^I"-1")-row.η))
end
function axial_constants(data,c,d,F)
    @assert sup(d)<data.profile.L-data.profile.d
    B,z2,z=F.B,F.z2,F.z
    kv=knot_interval(data.profile.tv,c)
    Bv=proof_float64(local_basis(data.profile.tv,data.profile.p,kv,c,d))
    return (;kv,Bv,hv=proof_float64((d-c)/I"2"),saxis=isthinzero(c),
             Zz1=reduce_order(I"6"*z2*B^I"-3"),Zzz1=reduce_order(I"36"*z2*B^I"-6"),
             Zzz2=reduce_order(I"9"*(B^4+I"1")*B^I"-7"),z=reduce_order(z))
end
function quick_field(data,R,Z,row,zeta)
    p=data.profile.p
    Cp=data.profile.Ψ[R.ku-p:R.ku,Z.kv-p:Z.kv]
    M=R.Bu'*Cp*Z.Bv
    P=zeros(Interval{Float64},Taylor(p)⊗Taylor(p))
    for i=0:p,j=0:p
        P[(i,j)]=M[i+1,j+1]
    end
    function deriv(i,j)
       Q=Model((one(first(coefficients(P)))/(R.hu^interval(i)*Z.hv^interval(j)))*differentiate(P,(i,j)),FIELD_ORDER)
       q=reduce_order(Q)
       return Z.saxis ? reduce_order(substitute_axis(q)) : q
    end
    p00,p10,p01,p11,p20,p02=deriv(0,0),deriv(1,0),deriv(0,1),deriv(1,1),deriv(2,0),deriv(0,2)
    psiz=p00+Z.Zz1*p01
    rpsirz=R.Rzr*(p10+Z.Zz1*p11)
    psirr=Z.z*(R.Rrr1*p20+R.Rrr2*p10)
    psizz=Z.z*(Z.Zzz1*p02+Z.Zzz2*p01)
    psir=R.Rr1*Z.z*p10
    return R.basea-rpsirz-I"2"*psiz,(R.r*(psirr-psizz)+I"3"*psir)/I"2",
        I"0.25"-I"0.100001"+rpsirz+psiz-row.η-zeta
end

function axial_moments(data,ctx,vertical,c,d)
    row=ctx.row
    modes=gram_modes(row); bases=data.coercivity.bases; bz=bases[row.z_basis]
    mz=sort!(unique(k for (_,k,_) in modes))
    zpairs=[(k,l) for k in mz for l in mz if k<=l]
    return get!(vertical,(inf(c),inf(d))) do
        a,b=ctx.U[1:2]
        F=coefficient_models(data.profile,(a,b,c,d),12,data.far)
        z=gram_order(F.z,24); Iz=strip_z_range(z,c,d)
        Z=Dict(k=>functional_mode(bz,k,z,Iz;N=24).exact.f for k in mz)
        Uz=functional_bump(z,Iz,bz.nodes[[1,3,4,6]]).f
        H=gram_moments(Uz*gram_order(F.dz,24),Z,2,FIELD_ORDER)
        return (;ax=axial_constants(data,c,d,F),H,Ht=transpose(column_table(H,zpairs,FIELD_ORDER)),h=[H.norms[p] for p in zpairs])
    end
end
function grams_child(data,row,a,b,V,vertical,zetas)
    modes=gram_modes(row);br=data.coercivity.bases[row.r_basis]
    mr=sort!(unique(m for (m,_,_) in modes)); mz=sort!(unique(k for (_,k,_) in modes))
    zpairs=[(k,l) for k in mz for l in mz if k<=l]
    ru=radial_constants(data,row,a,b,V); F1=ru.F
    r=gram_order(F1.r,24); T=functional_bump(r,F1.Ir,br.nodes[[1,3,4,6]]).f
    radial_basis=Dict(m=>functional_mode(br,m,r,F1.Ir;N=24).exact.f for m in mr)
    R=gram_moments(T*gram_order(F1.r2,24)*gram_order(F1.A,24)^3/I"6",radial_basis,1,FIELD_ORDER)
    S=[fill(I"0",FIELD_ORDER+1,length(zpairs)) for field in 1:6]; E=[fill(I"0",length(zpairs)) for field in 1:6]
    for j in eachindex(vertical)
        q=vertical[j]
        m,n,p=quick_field(data,ru,q.ax,row,zetas[j])
        jac=q.ax.saxis ? proof_float64((b-a)*sqrt(V[j+1])/I"4") : proof_float64((b-a)*(V[j+1]-V[j])/I"4")
        models=(m,n,p,m*m+n*n,n*(m+p),n*n+p*p)
        add_column_models!(S,E,models,q.H,zpairs,jac,FIELD_ORDER;Ht=q.Ht,h=q.h)
    end
    rp=column_pairs(R)
    return column_matrices(row,rp,zpairs,finish_column_models(R,rp,S,E,FIELD_ORDER))
end

# Production-order twelve Taylor contraction for the tightest enclosures.
function axial_moments12(data,ctx,vertical,c,d)
    row=ctx.row
    modes=gram_modes(row); bases=data.coercivity.bases; bz=bases[row.z_basis]
    mz=sort!(unique(k for (_,k,_) in modes))
    zpairs=[(k,l) for k in mz for l in mz if k<=l]
    return get!(vertical,(inf(c),inf(d))) do
        a,b=ctx.U[1:2]
        F=coefficient_models(data.profile,(a,b,c,d),12,data.far)
        z=gram_order(F.z,24); Iz=strip_z_range(z,c,d)
        Z=Dict(k=>functional_mode(bz,k,z,Iz;N=24).exact.f for k in mz)
        Uz=functional_bump(z,Iz,bz.nodes[[1,3,4,6]]).f
        H=gram_moments(Uz*gram_order(F.dz,24),Z,2,12)
        return (;z=gram_model_key(F.z),dz=gram_model_key(F.dz),H,Ht=transpose(column_table(H,zpairs,12)),h=[H.norms[p] for p in zpairs])
    end
end
function grams_child12(data,row,a,b,V,vertical,zetas)
    modes=gram_modes(row);br=data.coercivity.bases[row.r_basis]
    mr=sort!(unique(m for (m,_,_) in modes)); mz=sort!(unique(k for (_,k,_) in modes))
    zpairs=[(k,l) for k in mz for l in mz if k<=l]
    F1=coefficient_models(data.profile,(a,b,V[1],V[2]),12,data.far)
    radial_key=gram_model_key(F1.A)
    r=gram_order(F1.r,24); T=functional_bump(r,F1.Ir,br.nodes[[1,3,4,6]]).f
    radial_basis=Dict(m=>functional_mode(br,m,r,F1.Ir;N=24).exact.f for m in mr)
    R=gram_moments(T*gram_order(F1.r2,24)*gram_order(F1.A,24)^3/I"6",radial_basis,1,12)
    S=[fill(I"0",13,length(zpairs)) for field in 1:6]; E=[fill(I"0",length(zpairs)) for field in 1:6]
    for j in eachindex(vertical)
        F=j==1 ? F1 : coefficient_models(data.profile,(a,b,V[j],V[j+1]),12,data.far)
        q=vertical[j]
        gram_model_key(F.A)==radial_key && gram_model_key(F.z)==q.z && gram_model_key(F.dz)==q.dz || error("cached axial moments mismatch")
        m,n,p=matrix_field(F,row,zetas[j])
        models=(m,n,p,m*m+n*n,n*(m+p),n*n+p*p)
        add_column_models!(S,E,models,q.H,zpairs,F.jac,12;Ht=q.Ht,h=q.h)
    end
    rp=column_pairs(R)
    return column_matrices(row,rp,zpairs,finish_column_models(R,rp,S,E,12))
end

