# The integrals G,G₂ in (4.48). The coefficient models are two-dimensional;
# paired moments of the separated eigenfunctions avoid their full products.

function gram_integral(F::Model)
    all(isguaranteed,coefficients(F.Π)) && isguaranteed(F.ϱ) && inf(F.ϱ)≥0 || error("invalid integral model")
    S=I"0"
    for (i,j) in indices(space(F.Π))
        iseven(i) && iseven(j) || continue
        S=S+I"4"*F.Π[(i,j)]/interval(Float64,(i+1)*(j+1))
    end
    return S+I"[-4,4]"*F.ϱ
end

function gram_moment(F::Model,axis,α)
    S=I"0"
    for (i,j) in indices(space(F.Π))
        c=F.Π[(i,j)]
        isthinzero(c) && continue
        (axis==1 ? j==0 : i==0) || error("a separated moment depends on both coordinates")
        k=axis==1 ? i : j
        iseven(k+α) && (S=S+I"2"*c/interval(Float64,k+α+1))
    end
    return S+I"[-2,2]"*F.ϱ
end

# Raising the truncation order preserves the existing model and its error.
gram_order(F::Model,N)=Model(F.Π,F.ϱ,N)
gram_model_key(F::Model)=(order(space(F.Π)),Tuple((inf(c),sup(c)) for c in coefficients(F.Π)),inf(F.ϱ),sup(F.ϱ),F.N)

function gram_moments(weight,modes,axis,N)
    table=Dict{Tuple{Int,Int},Vector{Interval{Float64}}}()
    for m in sort(collect(keys(modes))),n in sort(collect(keys(modes)))
        m≤n || continue
        P=weight*modes[m]*modes[n]
        table[(m,n)]=[gram_moment(P,axis,α) for α in 0:N]
    end
    squares=Dict(m=>intersect_interval(table[(m,m)][1],I"[0,inf]") for m in keys(modes))
    all(x->!isempty_interval(x)&&isguaranteed(x),values(squares)) || error("a squared mode integral is negative")
    norms=Dict((m,n)=>sqrt(squares[m]*squares[n]) for (m,n) in keys(table))
    return (table=table,norms=norms)
end

function gram_pairing(F::Model,R,Z,rpair,zpair)
    a=R.table[rpair];b=Z.table[zpair];S=I"0"
    for (i,j) in indices(space(F.Π))
        S=S+F.Π[(i,j)]*a[i+1]*b[j+1]
    end
    return S+I"[-1,1]"*F.ϱ*R.norms[rpair]*Z.norms[zpair]
end

function gram_modes(row)
    return vcat([(m,k,1) for (m,k) in eachrow(row.sigma1)],
                [(m,k,2) for (m,k) in eachrow(row.sigma2)])
end

function gram_eigenvalue(b,m)
    return m==0 ? I"0" : proof_float64(b.modes[m].lambda)
end

# F.measure includes the normalized-cell Jacobian and the physical r³ measure.
# R,Z contain the corresponding one-dimensional factors, without F.jac.
function gram_cell(F,row,bases,ζ=I"0";N=24,radial=Dict(),vertical=Dict())
    modes=gram_modes(row);K=length(modes)
    K>0 || error("a retained cell space is empty")
    br=bases[row.r_basis];bz=bases[row.z_basis]
    r=gram_order(F.r,N);z=gram_order(F.z,N)
    Ir=F.Ir;Iz=functional_range(z)
    mr=sort(unique(m for (m,k,a) in modes));mz=sort(unique(k for (m,k,a) in modes))
    M=matrix_field(F,row,ζ)
    if row.s==1
        entries=(M[1],);squares=(M[1]*M[1],)
    else
        a,b,c=M
        entries=(a,b,c);squares=(a*a+b*b,b*(a+c),b*b+c*c)
    end
    order2=maximum(maximum(order(space(v.Π))) for v in (entries...,squares...))
    Rt=get!(radial,(row.r_basis,Tuple(mr),N,order2,gram_model_key(F.A))) do
        Tr=functional_bump(r,Ir,br.nodes[[1,3,4,6]]).f
        R=Dict(m=>functional_mode(br,m,r,Ir;N=N).exact.f for m in mr)
        ρr=Tr*gram_order(F.r2,N)*gram_order(F.A,N)^3/I"6"
        gram_moments(ρr,R,1,order2)
    end
    Zt=get!(vertical,(row.z_basis,Tuple(mz),N,order2,gram_model_key(F.B))) do
        Uz=functional_bump(z,Iz,bz.nodes[[1,3,4,6]]).f
        Z=Dict(k=>functional_mode(bz,k,z,Iz;N=N).exact.f for k in mz)
        gram_moments(Uz*gram_order(F.dz,N),Z,2,order2)
    end
    G=fill(I"0",K,K);G₂=copy(G)
    factor=(row.z_basis=="U1" ? I"1" : I"2")*F.jac
    for p in 1:K,q in p:K
        m,k,a=modes[p];n,l,b=modes[q]
        entry=row.s==1 ? 1 : a==b ? (a==1 ? 1 : 3) : 2
        rp=(min(m,n),max(m,n));zp=(min(k,l),max(k,l))
        G[p,q]=factor*gram_pairing(entries[entry],Rt,Zt,rp,zp)
        G₂[p,q]=factor*gram_pairing(squares[entry],Rt,Zt,rp,zp)
        G[q,p]=G[p,q];G₂[q,p]=G₂[p,q]
    end
    Λ=[gram_eigenvalue(br,m)+gram_eigenvalue(bz,k) for (m,k,a) in modes]
    all(A->all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),A),(G,G₂,Λ)) || error("an invalid Gram enclosure was produced")
    return (G=G,G2=G₂,Λ=Λ,modes=modes)
end

# The middle term of (4.62) is bounded by orthonormality and 0≤ρ≤1.
# Parity projection cancels κ² on a reflected pair.
function vector_B2(row)
    row.s==2 || error("a vector row is required")
    return [interval(Float64,f.nu)*sum((interval(Float64,c)^2 for c in (f.xi1...,f.xi2...));init=I"0") for f in row.fs]
end

# Interior of (4.63), reusing the exact radial primitive from Functionals.jl.
# Its normalization already includes ν, κ and the parity projection; 4π²
# restores the positive-half-domain physical integral.
function scalar_B2_cell(F,row,bases;N=24)
    row.s==1 || error("a scalar row is required")
    r=gram_order(F.r,N);z=gram_order(F.z,N)
    cache=Dict();B=Interval{Float64}[]
    for f in row.fs
        Q=functional_models(f,bases,r,z,F.Ir,functional_range(z);cache=cache,N=N)
        value=I"4"*interval(Float64,π)^2*gram_integral(F.measure*Q.Vr*Q.Vr)
        value=intersect_interval(value,I"[0,inf]")
        isempty_interval(value) && error("a squared scalar functional integral is negative")
        push!(B,value)
    end
    return B
end

# H(rb,z) contains only the constant radial mode: the other modes have zero
# weighted mean. Since U²≤U, vertical orthonormality bounds the infinite tail.
function scalar_B2_tail(f,bases)
    f.kind=="scalar" || error("a scalar functional is required")
    br=bases[f.r_basis];rb=proof_float64(br.nodes[end])
    X=Model(Sequence(Taylor(0)⊗Taylor(0),[I"1"]),24)
    mass=functional_bump_integral(X,f.r_spec,4).Π[(0,0)]
    S=I"0"
    for n in eachindex(f.xi1)
        f.sigma1[n,1]==0 || continue
        c=interval(Float64,f.xi1[n]);S=S+c*c
    end
    return interval(Float64,f.nu)*mass*S/(I"2"*rb*rb)
end
