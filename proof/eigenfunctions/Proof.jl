# The same recurrence, matching, zero-count and normalization estimates as the
# existing Appendix A proof. The normalized target is fixed by the shared data.
function weight_context(name,roots)
    _,strings,radial,M=only(filter(w->w[1]==name,WEIGHTS))
    a,x2,x3,b=parse.(Interval{BigFloat},strings)
    x1=(a+x2)/interval(2);x4=(x3+b)/interval(2)
    nodes=(a,x1,x2,x3,x4,b)
    pieces=series_pieces(nodes...,radial,roots,ν,N)
    h=[x1-a,(x2-x1)/interval(2),(x4-x3)/interval(2),b-x4]
    c=[interval(BigFloat,1),inv(exp(interval(BigFloat,1))*log(ν)),
       (interval(2)/(exp(interval(BigFloat,1))*log(ν)))^2]
    return (name=name,radial=radial,nodes=nodes,pieces=pieces,h=h,c=c,M=M)
end

function validate_approximation(r,name,m)
    r["N"]==N && r["prec"]==128 && r["m"]==m || error("incorrect eigenpair input: $name $m")
    length(r["theta"])==9 || error("the matching vector must have nine entries")
    for a in (r["theta"],coefficients(r["phi_a"]),coefficients(r["psi_a"]),coefficients(r["psi_b"]),coefficients(r["phi_b"]))
        all(x->x isa Interval{BigFloat}&&isguaranteed(x)&&isthin(x)&&isfinite(inf(x))&&isfinite(sup(x))&&precision(inf(x))==precision(sup(x))==128,a) || error("input coefficients must be finite thin 128-bit intervals")
    end
    all(length(coefficients(r[k]))==N+1 for k in ("phi_a","psi_a","psi_b","phi_b")) || error("incorrect series lengths")
    return true
end

function prove_mode(context,m,fixed)
    name=context.name;radial=context.radial
    a,x1,x2,x3,x4,b=context.nodes
    φa,ψa,ψb,φb=context.pieces;c=context.c;h=context.h
    input=load(joinpath(DATA,"Eigenfunctions",mode_file(name,m)))
    validate_approximation(input,name,m)
    θ̄=input["theta"];λ̄=θ̄[1]
    isequal_interval(λ̄,fixed.lambda_bar) || error("the delivered approximation uses another eigenvalue")
    φ̃=(input["phi_a"],input["psi_a"],input["psi_b"],input["phi_b"])
    ψ̃=Tuple(Sequence(Taylor(N),getproperty(fixed,key)) for key in (:phi_a,:psi_a,:psi_b,:phi_b))
    all(x->isguaranteed(x)&&isthin(x)&&isfinite(inf(x)),vcat([coefficients(p) for p in ψ̃]...,[fixed.A,fixed.B])) || error("invalid normalized target")
    e9=[fill(interval(BigFloat,0),8);interval(BigFloat,1)]

    ops=solution_operators(φa,ψa,ψb,φb,λ̄,ν,N)
    L,dL=computation_L(φa,ψa,ψb,φb,ops,radial,x2,x3,λ̄,ν)
    M=interval.(inv(mid.([dL*θ̄[2:9] L])))
    ε=maximum(abs.(M*(L*θ̄[2:9]-e9)))
    Λ=λ̄+interval(-sup(rstar),sup(rstar));θ̃=θ̄[2:9].+interval(-sup(rstar),sup(rstar))
    ops=solution_operators(φa,ψa,ψb,φb,Λ,ν,N)
    L,dL=computation_L(φa,ψa,ψb,φb,ops,radial,x2,x3,Λ,ν)
    Z=maximum(sum(abs.(UniformScaling(interval(BigFloat,1))-M*[dL*θ̃ L]),dims=2))
    isguaranteed(ε)&&isguaranteed(Z)&&sup(Z)<1 || error("the matching contraction failed: $name $m")
    r=ε/(interval(1)-Z)
    sup(r)<=inf(rstar) || error("the matching radius failed: $name $m")
    Λr=λ̄+interval(-sup(r),sup(r));Θ=θ̄.+interval(-sup(r),sup(r))

    Sa,Ca10,Ca01,Cb10,Cb01,Sb=solution_operators(φa,ψa,ψb,φb,Λr,ν,N)
    e=[norm(Θ[2]*Sa[1]-φ̃[1],X)+abs(Θ[2])*Sa[2],
       norm(Θ[3]*Ca10[1]+Θ[4]*Ca01[1]-φ̃[2],X)+abs(Θ[3])*Ca10[2]+abs(Θ[4])*Ca01[2],
       norm(Θ[7]*Cb10[1]+Θ[8]*Cb01[1]-φ̃[3],X)+abs(Θ[7])*Cb10[2]+abs(Θ[8])*Cb01[2],
       norm(Θ[9]*Sb[1]-φ̃[4],X)+abs(Θ[9])*Sb[2]]
    S_i,D_i=plateau_bounds(radial,x2,x3,Λr)
    Δ=[max(maximum([c[i+1]*e[j]/pown(h[j],i) for j=1:4]),
           r*(interval(2)*S_i[i+1]+(abs(θ̄[5])+abs(θ̄[6]))*D_i[i+1])) for i=0:2]
    all(isguaranteed,Δ) || error("unguaranteed derivative errors")
    nz,zeros_in,supφ̃=count_zeros(context.nodes,radial,context.pieces,φ̃,λ̄,θ̄[5],θ̄[6],sup(max(Δ[1],Δ[2])))
    nz==m || error("the zero count is $nz instead of $m for $name")
    𝒩=computation_𝒩(context.pieces,φ̃,e,supφ̃,radial,x2,x3,Λr,Θ[5],Θ[6])
    isguaranteed(𝒩)&&inf(𝒩)>0 || error("positive normalization failed")
    n𝒩=inv(sqrt(𝒩))
    normalized=Tuple(n𝒩*p for p in φ̃)
    A=n𝒩*θ̄[5];B=n𝒩*θ̄[6]
    if isequal_interval(a,-b)&&isequal_interval(x2,-x3)
        normalized,A,B=symmetrise(normalized,A,B,(-1)^m)
    end
    d=[norm(normalized[j]-ψ̃[j],X) for j=1:4]
    Δn=[abs(n𝒩)*Δ[i+1]+max(maximum([c[i+1]*d[j]/pown(h[j],i) for j=1:4]),
           (abs(A-fixed.A)+abs(B-fixed.B))*S_i[i+1]) for i=0:2]
    if isequal_interval(a,-b)&&isequal_interval(x2,-x3)
        # The delivered polynomial and plateau target must already have parity (-1)^m.
        projected,A,B=symmetrise(ψ̃,fixed.A,fixed.B,(-1)^m)
        all(all(isequal_interval.(coefficients(p),coefficients(q))) for (p,q) in zip(projected,ψ̃)) &&
            isequal_interval(A,fixed.A)&&isequal_interval(B,fixed.B) || error("the fixed target has incorrect parity")
    end
    result=(name=name,m=m,N=N,prec=128,radial=radial,nodes=collect(context.nodes),
            lambda=outward128(Λr),lambda_bar=λ̄,epsilon=outward128(ε),Z=outward128(Z),r=outward128(r),theta=outward128.(Θ),
            phi_a=fixed.phi_a,psi_a=fixed.psi_a,psi_b=fixed.psi_b,phi_b=fixed.phi_b,
            A=fixed.A,B=fixed.B,errors=outward128.(Δn),zeros=nz,
            zero_intervals=[interval(BigFloat,lo,hi) for (lo,hi) in zeros_in],norm2=outward128(𝒩))
    check_mode(result,name,m)
    return result
end

# Exact weighted mass of the plateau bump, including both polynomial ramps.
function constant_mode(name)
    _,strings,radial,M=only(filter(w->w[1]==name,WEIGHTS))
    setprecision(512) do
        a,x2,x3,b=parse.(Interval{BigFloat},strings)
        if radial
            ramp(c,h)=h*sum(interval(BigFloat,binomial(3,k))*pown(c,3-k)*pown(h,k)*
                         (interval(BigFloat,126)/interval(BigFloat,k+6)-interval(BigFloat,420)/interval(BigFloat,k+7)+
                          interval(BigFloat,540)/interval(BigFloat,k+8)-interval(BigFloat,315)/interval(BigFloat,k+9)+interval(BigFloat,70)/interval(BigFloat,k+10)) for k=0:3)
            mass=(x3^4-x2^4)/interval(BigFloat,4)+ramp(a,x2-a)-ramp(b,x3-b)
        else
            mass=(x3-x2)+(x2-a+b-x3)/interval(BigFloat,2)
        end
        inf(mass)>0 || error("the constant mode has nonpositive mass")
        # Keep a guaranteed interval; conversion rounds outward to 128 bits.
        value=inv(sqrt(mass))
        return setprecision(128) do
            outward128(value)
        end
    end
end
