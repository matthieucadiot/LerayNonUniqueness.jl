# Integrate one radial column before contracting the radial paired moments.
# For each matrix entry, S=Σ_j jac_j P_j H_jᵀ. The radial table multiplies
# S only once, after all vertical cells; the same applies to Cauchy errors.

# Keep the output order while evaluating independent column states in parallel.
function threaded_map(f,items)
    value=f(first(items))
    results=Vector{typeof(value)}(undef,length(items));results[1]=value
    Threads.@threads for k in 2:length(items)
        results[k]=f(items[k])
    end
    return results
end

function column_coefficients(F::Model,N)
    P=fill(I"0",N+1,N+1)
    for (i,j) in indices(space(F.Π))
        i≤N && j≤N || error("a coefficient model exceeds the column order")
        P[i+1,j+1]=F.Π[(i,j)]
    end
    return P
end

column_pairs(T)=sort(collect(keys(T.table)))
column_table(T,pairs,N)=[T.table[p][j+1] for p in pairs,j in 0:N]

function add_column_models!(S,E,models,H,pairs,jac,N;Ht=transpose(column_table(H,pairs,N)),h=[H.norms[p] for p in pairs])
    isguaranteed(jac) && inf(jac)>0 || error("the integration Jacobian must be positive and guaranteed")
    for f in eachindex(models)
        M=models[f]
        isguaranteed(M.ϱ) && inf(M.ϱ)≥0 && isfinite(sup(M.ϱ)) || error("invalid coefficient remainder")
        # Exact spline terms can have degree above the product truncation order.
        # Charge every discarded coefficient to the uniform remainder first.
        M=functional_order(M,N)
        S[f].+=jac.*(column_coefficients(M,N)*Ht)
        E[f].+=(jac*M.ϱ).*h
    end
    return nothing
end

function finish_column_models(R,pairs,S,E,N)
    Rt=column_table(R,pairs,N);r=[R.norms[p] for p in pairs]
    return [Rt*S[f]+[I"[-1,1]"*r[p]*E[f][q] for p in eachindex(r),q in eachindex(E[f])] for f in eachindex(S)]
end

function column_matrices(row,rpairs,zpairs,integrals)
    modes=gram_modes(row);K=length(modes)
    ir=Dict(p=>i for (i,p) in enumerate(rpairs))
    iz=Dict(p=>i for (i,p) in enumerate(zpairs))
    G=fill(I"0",K,K);G2=copy(G)
    factor=row.z_basis=="U1" ? I"1" : I"2"
    offset=row.s==1 ? 1 : 3
    for p in 1:K,q in p:K
        m,k,a=modes[p];n,l,b=modes[q]
        f=row.s==1 ? 1 : a==b ? (a==1 ? 1 : 3) : 2
        r=ir[(min(m,n),max(m,n))];z=iz[(min(k,l),max(k,l))]
        G[p,q]=factor*integrals[f][r,z]
        G2[p,q]=factor*integrals[f+offset][r,z]
        G[q,p]=G[p,q];G2[q,p]=G2[p,q]
    end
    return G,G2
end

# P_m=H_m/r³ in (4.63). The nonconstant modes use their exact differential
# equation; the constant mode uses the polynomial primitive of t³T(t).
function scalar_radial_matrix(F,br,R,mr,N)
    r=gram_order(F.r,N);Ir=F.Ir;spec=br.nodes[[1,3,4,6]]
    T=functional_bump(r,Ir,spec).f
    P=Dict(m=>(m==0 ? proof_float64(br.constant)*(functional_radial_primitive(r,Ir,spec)*r^I"-3") :
               -(I"1"/proof_float64(br.modes[m].lambda))*(T*R[m].d)) for m in mr)
    weight=gram_order(F.r2,N)*gram_order(F.A,N)^3/I"6"
    J=gram_moments(weight,P,1,0)
    return [J.table[(min(m,n),max(m,n))][1] for m in mr,n in mr]
end

function scalar_column_bounds(row,mr,mz,J,K)
    ir=Dict(m=>i for (i,m) in enumerate(mr));iz=Dict(k=>i for (i,k) in enumerate(mz))
    B=Interval{Float64}[]
    for f in row.fs
        C=fill(I"0",length(mr),length(mz))
        for n in eachindex(f.xi1)
            m,k=f.sigma1[n,:]
            C[ir[m],iz[k]]+=interval(Float64,f.xi1[n])
        end
        value=interval(Float64,f.nu)*sum(C.*(J*C*K))
        value=intersect_interval(value,I"[0,inf]")
        isempty_interval(value) && error("a scalar column integral is negative")
        push!(B,value)
    end
    return B
end

# G,G2 and scalar B² are contributions of the indicated radial columns.
# Keep one set of vertical moments in memory and reuse it across the batch.
# The assembler uses vector_B2 once, and adds scalar_B2_tail once.
function gram_columns(data,row,bases,B,far,U,V,columns;N=24,progress=false,on_column=nothing)
    columns=collect(columns)
    !isempty(columns) && length(unique(columns))==length(columns) &&
        all(i->1≤i<length(U),columns) && length(V)≥2 && N≥12 || error("invalid columns or Taylor order")
    N2=12;modes=gram_modes(row)
    br=bases[row.r_basis];bz=bases[row.z_basis]
    mr=sort(unique(m for (m,k,a) in modes));mz=sort(unique(k for (m,k,a) in modes))
    zpairs=[(k,l) for k in mz for l in mz if k≤l]
    nf=row.s==1 ? 2 : 6
    scalar=row.s==1 && row.n>0
    states=threaded_map(columns) do i
        F=coefficient_models(data,(U[i],U[i+1],V[1],V[2]),N2,far)
        r=gram_order(F.r,N);Ir=F.Ir
        R=Dict(m=>functional_mode(br,m,r,Ir;N=N).exact for m in mr)
        T=functional_bump(r,Ir,br.nodes[[1,3,4,6]]).f
        weight=T*gram_order(F.r2,N)*gram_order(F.A,N)^3/I"6"
        Rt=gram_moments(weight,Dict(m=>R[m].f for m in mr),1,N2)
        S=[fill(I"0",N2+1,length(zpairs)) for f in 1:nf]
        E=[fill(I"0",length(zpairs)) for f in 1:nf]
        J=scalar ? scalar_radial_matrix(F,br,R,mr,N) : nothing
        K=fill(I"0",length(mz),length(mz))
        return (i=i,F=F,radial_key=gram_model_key(F.A),Rt=Rt,S=S,E=E,J=J,K=K)
    end
    for j in 1:length(V)-1
        first_state=first(states);i=first_state.i
        first_F=j==1 ? first_state.F : coefficient_models(data,(U[i],U[i+1],V[j],V[j+1]),N2,far)
        vertical_key=gram_model_key(first_F.B)
        z=gram_order(first_F.z,N);Iz=functional_range(z)
        Z=Dict(k=>functional_mode(bz,k,z,Iz;N=N).exact.f for k in mz)
        Uz=functional_bump(z,Iz,bz.nodes[[1,3,4,6]]).f
        dz=gram_order(first_F.dz,N)
        H=gram_moments(Uz*dz,Z,2,N2)
        H2=scalar ? gram_moments(Uz*Uz*dz,Z,2,0) : nothing
        Ht=transpose(column_table(H,zpairs,N2));h=[H.norms[p] for p in zpairs]
        Threads.@threads for index in eachindex(states)
            let state=states[index], i=states[index].i
            box=(U[i],U[i+1],V[j],V[j+1])
            F=i==first_state.i ? first_F : j==1 ? state.F : coefficient_models(data,box,N2,far)
            gram_model_key(F.A)==state.radial_key || error("the radial model changes inside its column")
            gram_model_key(F.B)==vertical_key || error("the vertical model changes inside its row")
            ζ=row.block==:odd ? zeta_bound(B,box) : I"0"
            M=matrix_field(F,row,ζ)
            models=if row.s==1
                (M[1],M[1]*M[1])
            else
                a,b,c=M
                (a,b,c,a*a+b*b,b*(a+c),b*b+c*c)
            end
            add_column_models!(state.S,state.E,models,H,zpairs,F.jac,N2;Ht=Ht,h=h)
            if scalar
                # K includes the full cell Jacobian, including its radial factor.
                # Consequently J has no additional radial Jacobian.
                for (p,k) in enumerate(mz),(q,l) in enumerate(mz)
                    state.K[p,q]+=I"2"*F.jac*H2.table[(min(k,l),max(k,l))][1]
                end
            end
            end
        end
        if progress && (j%25==0 || j==length(V)-1)
            println(row.id,": vertical cell ",j,"/",length(V)-1,"; columns ",first(columns),"..",last(columns))
            flush(stdout)
        end
    end
    Λ=[gram_eigenvalue(br,m)+gram_eigenvalue(bz,k) for (m,k,a) in modes]
    callback_lock=ReentrantLock()
    return threaded_map(states) do state
        rpairs=column_pairs(state.Rt)
        integrals=finish_column_models(state.Rt,rpairs,state.S,state.E,N2)
        G,G2=column_matrices(row,rpairs,zpairs,integrals)
        B2=scalar ? scalar_column_bounds(row,mr,mz,state.J,state.K) : fill(I"0",row.n)
        all(A->all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),A),(G,G2,Λ,B2)) || error("an invalid column integral was produced")
        result=(column=state.i,G=G,G2=G2,Λ=Λ,B2=B2,cells=length(V)-1)
        if on_column!==nothing
            lock(callback_lock) do
                on_column(result)
            end
        end
        return result
    end
end

function gram_column(data,row,bases,B,far,U,V,i;N=24,progress=false)
    return first(gram_columns(data,row,bases,B,far,U,V,(i,);N=N,progress=progress))
end
