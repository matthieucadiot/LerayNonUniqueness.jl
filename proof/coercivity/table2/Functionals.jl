# Functionals.jl: the physical q_j of (4.64), their numerical versions, and the vector fields
# used in (4.55), (4.63). The existing Model arithmetic carries all coefficient and function errors.

function functional_constant(c,X)
    P=zeros(Interval{Float64},Taylor(0) ⊗ Taylor(0)); P[(0,0)]=interval(Float64,c)
    return Model(P,X.N)
end

function functional_range(M)
    e=M.ϱ
    for (i,j) in indices(space(M.Π))
        (i==0 && j==0) || (e=e+abs(M.Π[(i,j)]))
    end
    return M.Π[(0,0)]+I"[-1,1]"*e
end

# A varying function in an interval is represented by a thin constant plus a uniform error.
# A wide constant coefficient alone would not justify integration as a polynomial.
function functional_range_model(v,X)
    c=interval(Float64,Float64(mid(v)))
    return Model(functional_constant(c,X).Π,interval(Float64,sup(abs(v-c))),X.N)
end

function functional_order(X,N)
    P,e=πN(X.Π,min(X.N,N))
    return Model(P,X.ϱ+e,min(X.N,N))
end

function functional_round(x,prec)
    isguaranteed(x) || error("a functional coefficient is not guaranteed")
    return interval(Float64,x)
end

# DLMF 10.2.2 and 10.8.1, n=0,1. The omitted ordinary and harmonic-weighted
# series are bounded by geometric sums, after their absolute term ratio is <1.
function functional_bessel(n,x)
    n in (0,1) && inf(x)>0 || error("Bessel seeds require n=0 or1 and x>0")
    precision(BigFloat)==1024 || error("Bessel seeds require 1024-bit working precision")
    w=x*x/I"4"; t=n==0 ? one(x) : x/I"2"
    H=zero(x); γ=RESPONSE_BESSEL_EULER
    SJ=zero(x); SY=zero(x)
    tolerance=big(I"2")^I"-700"
    for k in 0:4095
        kn=interval(BigFloat,k+1)
        h=I"2"*H+(n==1 ? one(x)/kn : zero(x))-I"2"*γ
        SJ=SJ+t; SY=SY+h*t
        tn=-t*w/(kn*interval(BigFloat,k+n+1))
        Hn=H+one(x)/kn
        hn=I"2"*Hn+(n==1 ? one(x)/interval(BigFloat,k+2) : zero(x))-I"2"*γ
        q=abs(w)/(interval(BigFloat,k+2)*interval(BigFloat,k+n+2))
        if sup(q)<1
            TJ=abs(tn)/(I"1"-q)
            TY=abs(tn)*(abs(hn)/(I"1"-q)+I"2"*q/(interval(BigFloat,k+2)*(I"1"-q)^2))
            if sup(max(TJ,TY))<inf(tolerance)
                J=SJ+I"[-1,1]"*TJ
                πb=RESPONSE_BESSEL_PI
                finite=n==0 ? zero(x) : -I"2"/(πb*x)
                Y=finite+I"2"/πb*log(x/I"2")*J-(SY+I"[-1,1]"*TY)/πb
                return J,Y
            end
        end
        t=tn; H=Hn
    end
    error("Bessel seed series did not reach its certified tail bound")
end

const functional_seed_cache=Dict{Tuple,Tuple}()
const functional_seed_lock=ReentrantLock()

function functional_central_seed(b,m,c)
    key=(b.radial,inf(m.lambda_bar),sup(m.lambda_bar),inf(m.A),sup(m.A),inf(m.B),sup(m.B),inf(c),sup(c))
    saved=lock(functional_seed_lock) do
        get(functional_seed_cache,key,nothing)
    end
    saved===nothing || return saved
    prec=precision(BigFloat)
    a0,a1=setprecision(1024) do
        unit=big(I"1"); λ=unit*m.lambda_bar; cc=unit*c
        k=sqrt(λ); A=unit*m.A; B=unit*m.B
        if b.radial
            J0,Y0=functional_bessel(0,k*cc)
            J1,Y1=functional_bessel(1,k*cc)
            g=(A*J1+B*Y1)/cc
            gd=k*(A*J0+B*Y0)/cc-I"2"*g/cc
            return g,gd
        else
            return A*cos(k*cc)+B*sin(k*cc),k*(-A*sin(k*cc)+B*cos(k*cc))
        end
    end
    value=(functional_round(a0,prec),functional_round(a1,prec))
    return lock(functional_seed_lock) do
        get!(functional_seed_cache,key,value)
    end
end

# On the central plateau the numerical mode solves g''+3g'/r+lambda_bar*g=0
# (radial) or g''+lambda_bar*g=0 (axial). Its Taylor coefficients are obtained
# algebraically. A complex-disc Gronwall bound and Cauchy's estimate bound the tail.
function functional_central(b,m,X,I)
    c=interval(Float64,Float64(mid(I))); d=interval(Float64,sup(abs(I-c)))
    a0,a1=functional_central_seed(b,m,c)
    λ=interval(Float64,m.lambda_bar); s=max(I"1",sqrt(λ))
    R=b.radial ? min(c/I"2",interval(Float64,X.N+1)/s) : max(I"2"*d,interval(Float64,X.N+1)/s)
    R=interval(Float64,inf(R))
    inf(R)>0 && sup(d)<inf(R) || error("subdivide the functional cell: the Taylor disc is too small")
    M=b.radial ? max(s,λ/s+I"3"/(c-R)) : max(s,λ/s)
    B=max(abs(a0),abs(a1)/s)*exp(M*R)
    a=[a0,a1]
    for n in 0:X.N
        den=interval(Float64,(n+2)*(n+1))
        if b.radial
            anm1=n==0 ? zero(a0) : a[n]
            push!(a,-(interval(Float64,(n+1)*(n+3))*a[n+2]+λ*c*a[n+1]+λ*anm1)/(c*den))
        else
            push!(a,-λ*a[n+1]/den)
        end
    end
    q=d/R
    tail0=B*q^interval(Float64,X.N+1)/(I"1"-q)
    tail1=B/R*interval(Float64,X.N+2)*q^interval(Float64,X.N+1)/(I"1"-q)^2
    tail2=B/(R*R)*interval(Float64,(X.N+3)*(X.N+2))*q^interval(Float64,X.N+1)*(I"1"+q)/(I"1"-q)^3
    Y=X-c
    jets=Model[]
    for derivative in 0:2
        coefficients=[a[j+derivative+1]*interval(Float64,derivative==0 ? 1 : derivative==1 ? j+1 : (j+1)*(j+2)) for j in 0:X.N]
        P=functional_constant(coefficients[end],X)
        for j in X.N-1:-1:0
            P=coefficients[j+1]+Y*P
        end
        ε=(tail0,tail1,tail2)[derivative+1]
        push!(jets,Model(P.Π,P.ϱ+ε,P.N))
    end
    return (f=jets[1],d=jets[2],dd=jets[3])
end

function functional_polynomial_jet(coefficients,Y,scale)
    n=length(coefficients)-1
    jets=Model[]
    for derivative in 0:2
        if n<derivative
            push!(jets,functional_constant(I"0",Y)); continue
        end
        c=[interval(Float64,coefficients[j+derivative+1])*interval(Float64,derivative==0 ? 1 : derivative==1 ? j+1 : (j+1)*(j+2)) for j in 0:n-derivative]
        P=functional_constant(c[end],Y)
        for j in length(c)-1:-1:1
            P=c[j]+Y*P
        end
        derivative==1 && (P=scale*P)
        derivative==2 && (P=scale*scale*P)
        push!(jets,P)
    end
    return (f=jets[1],d=jets[2],dd=jets[3])
end

function functional_mode_piece(b,m,X,I,piece)
    a,x1,x2,x3,x4,bb=b.nodes
    if piece==1
        h=x1-a
        return functional_polynomial_jet(m.phi_a,(X-a)/h,I"1"/h)
    elseif piece==2
        c=(x1+x2)/I"2";h=(x2-x1)/I"2"
        return functional_polynomial_jet(m.psi_a,(X-c)/h,I"1"/h)
    elseif piece==3
        return functional_central(b,m,X,I)
    elseif piece==4
        c=(x3+x4)/I"2";h=(x4-x3)/I"2"
        return functional_polynomial_jet(m.psi_b,(X-c)/h,I"1"/h)
    else
        h=bb-x4
        return functional_polynomial_jet(m.phi_b,(bb-X)/h,-I"1"/h)
    end
end

function functional_mode(b,mindex,X,I;N=24)
    X=functional_order(X,N)
    domain=interval(inf(b.nodes[1]),sup(b.nodes[end]))
    J=intersect_interval(I,domain)
    isempty_interval(J) && error("functional mode requested outside its support")
    if inf(I)<inf(domain) || sup(I)>sup(domain)
        # Only the portion inside the plateau support is used; the multiplying bump
        # is zero elsewhere. Keep a range enclosure rather than extending a polynomial patch.
        v=functional_mode(b,mindex,functional_range_model(J,X),J;N=N)
        return (num=(; (k=>functional_range_model(functional_range(getproperty(v.num,k)),X) for k in (:f,:d,:dd))...),
                exact=(; (k=>functional_range_model(functional_range(getproperty(v.exact,k)),X) for k in (:f,:d,:dd))...))
    end
    if mindex==0
        z=functional_constant(I"0",X)
        v=functional_constant(b.constant,X)
        return (num=(f=v,d=z,dd=z),exact=(f=v,d=z,dd=z))
    end
    1≤mindex≤length(b.modes) || error("functional mode index $mindex is unavailable")
    m=b.modes[mindex]
    m.m==mindex || error("mode numbering is inconsistent")
    isguaranteed(m.lambda) && isguaranteed(m.lambda_bar) && inf(m.lambda)>0 && inf(m.lambda_bar)>0 || error("a positive certified mode eigenvalue is required")
    length(m.errors)==3 && all(x->isfinite(x)&&x≥0,m.errors) || error("invalid spectral derivative-error bounds")
    inf(I)≥inf(b.nodes[1]) && sup(I)≤sup(b.nodes[end]) || error("functional mode requested outside its support")
    candidates=Tuple{Int,Interval}[]
    for p in 1:5
        lo=max(inf(I),inf(b.nodes[p]));hi=min(sup(I),sup(b.nodes[p+1]))
        lo≤hi && push!(candidates,(p,interval(lo,hi)))
    end
    isempty(candidates) && error("functional cell has no spectral piece")
    if length(candidates)==1
        num=functional_mode_piece(b,m,X,I,candidates[1][1])
    else
        values=[functional_mode_piece(b,m,functional_range_model(J,X),J,p) for (p,J) in candidates]
        num=(; (k=>functional_range_model(reduce(hull,(functional_range(getproperty(v,k)) for v in values)),X) for k in (:f,:d,:dd))...)
    end
    exact=(f=Model(num.f.Π,num.f.ϱ+interval(Float64,m.errors[1]),num.f.N),
           d=Model(num.d.Π,num.d.ϱ+interval(Float64,m.errors[2]),num.d.N),
           dd=Model(num.dd.Π,num.dd.ϱ+interval(Float64,m.errors[3]),num.dd.N))
    return (num=num,exact=exact)
end

# A second, analytic upper bound from orthonormality and 0<=rho<=1. Integration
# of the vector models is preferred; this bound is useful as an independent check.
function functional_vector_bound(f)
    f.kind=="vector" || error("the orthonormal vector bound does not apply to scalar q")
    s=I"0"
    for c in (f.xi1...,f.xi2...)
        cc=interval(Float64,c);s=s+cc*cc
    end
    projection=f.symmetric_cell ? I"1" : I"0.5"
    return sqrt(interval(Float64,f.kappa)^2*interval(Float64,f.nu)*projection*s)
end

function validate_functional_parity(f)
    f.parity in (-1,1) || error("functional parity must be -1 or1")
    if f.symmetric_cell
        for (sigma,xi,target) in ((f.sigma1,f.xi1,f.parity),(f.sigma2,f.xi2,-f.parity))
            for n in eachindex(xi)
                (inf(xi[n])==0 && sup(xi[n])==0) && continue
                (iseven(sigma[n,2]) ? 1 : -1)==target || error("functional modes violate the prescribed component parity")
            end
        end
    end
    return true
end

# P[a,n1,n2,b] and its first two physical derivatives, Appendix A.
function functional_bump_piece(X,spec,piece)
    a,n1,n2,b=spec
    z=functional_constant(I"0",X)
    if piece==0
        return (f=z,d=z,dd=z)
    elseif piece==2
        return (f=functional_constant(I"1",X),d=z,dd=z)
    end
    h=piece==1 ? n1-a : b-n2
    y=piece==1 ? (X-a)/h : (b-X)/h
    scale=piece==1 ? I"1"/h : -I"1"/h
    c=[I"0",I"0",I"0",I"0",I"0",I"126",-I"420",I"540",-I"315",I"70"]
    return functional_polynomial_jet(c,y,scale)
end

function functional_bump(X,I,spec)
    a,n1,n2,b=spec
    if sup(I)≤inf(a) || inf(I)≥sup(b)
        return functional_bump_piece(X,spec,0)
    elseif inf(I)≥sup(n1) && sup(I)≤inf(n2)
        return functional_bump_piece(X,spec,2)
    elseif inf(I)≥sup(a) && sup(I)≤inf(n1)
        return functional_bump_piece(X,spec,1)
    elseif inf(I)≥sup(n2) && sup(I)≤inf(b)
        return functional_bump_piece(X,spec,3)
    end
    bounds=[functional_bump_piece(X,spec,0)]
    for (piece,lo,hi) in ((1,inf(a),sup(n1)),(2,inf(n1),sup(n2)),(3,inf(n2),sup(b)))
        l=max(inf(I),lo);h=min(sup(I),hi)
        l≤h && push!(bounds,functional_bump_piece(functional_range_model(interval(l,h),X),spec,piece))
    end
    return (; (k=>functional_range_model(reduce(hull,(functional_range(getproperty(v,k)) for v in bounds)),X) for k in (:f,:d,:dd))...)
end

# Integral of t^3 times the radial plateau bump. The polynomial is integrated
# after an affine substitution, with all binomial factors enclosed before arithmetic.
function functional_bump_integral(X,spec,piece)
    a,n1,n2,b=spec
    function ramp(Y,origin,h)
        total=functional_constant(I"0",Y)
        for (j,c) in ((5,I"126"),(6,-I"420"),(7,I"540"),(8,-I"315"),(9,I"70")), k in 0:3
            scalar=interval(Float64,binomial(3,k))*origin^interval(Float64,3-k)*h^interval(Float64,k+1)*c/interval(Float64,j+k+1)
            total=total+scalar*Y^(j+k+1)
        end
        return total
    end
    h1=n1-a;h2=b-n2
    left=ramp(functional_constant(I"1",X),a,h1)
    plateau=(n2^I"4"-n1^I"4")/I"4"
    right=-ramp(functional_constant(I"1",X),b,-h2)
    if piece==0
        return functional_constant(I"0",X)
    elseif piece==1
        return ramp((X-a)/h1,a,h1)
    elseif piece==2
        return left+(X^4-n1^I"4")/I"4"
    elseif piece==3
        return left+plateau+right+ramp((b-X)/h2,b,-h2)
    else
        return left+plateau+right
    end
end

function functional_radial_primitive(X,I,spec)
    a,n1,n2,b=spec
    for (piece,lo,hi) in ((0,-Inf,inf(a)),(1,sup(a),inf(n1)),(2,sup(n1),inf(n2)),(3,sup(n2),inf(b)),(4,sup(b),Inf))
        inf(I)≥lo && sup(I)≤hi && return functional_bump_integral(X,spec,piece)
    end
    values=Model[]
    for (piece,lo,hi) in ((0,-Inf,sup(a)),(1,inf(a),sup(n1)),(2,inf(n1),sup(n2)),(3,inf(n2),sup(b)),(4,inf(b),Inf))
        l=max(inf(I),lo);h=min(sup(I),hi)
        l≤h && push!(values,functional_bump_integral(functional_range_model(interval(l,h),X),spec,piece))
    end
    return functional_range_model(reduce(hull,(functional_range(v) for v in values)),X)
end

function functional_models(f,bases,r,z,Ir,Iz;cache=Dict(),N=24)
    validate_functional_parity(f)
    r=functional_order(r,N);z=functional_order(z,N)
    zeroM=functional_constant(I"0",r)
    # A rounded support boundary can leave a tiny overlap in the physical r=0
    # cell. The formula q=-div(V) contains1/r, but its compact support stays away
    # from the axis. Evaluate only that positive support portion and hull with0
    # for the remaining points. This also applies to the scalar primitive below a.
    if inf(Ir)≤0 && sup(Ir)>inf(f.r_spec[1])
        J=interval(inf(f.r_spec[1]),sup(Ir))
        inf(J)>0 || error("functional radial supports must avoid the axis")
        v=functional_models(f,bases,functional_range_model(J,r),z,J,Iz;cache=Dict(),N=N)
        fields=(; (k=>functional_range_model(hull(I"0",functional_range(getproperty(v,k))),r) for k in (:q,:qnum,:Vr,:Vz))...)
        difference=fields.q-fields.qnum
        return (q=fields.q,qnum=fields.qnum,δq=Model(zeroM.Π,sup_bound(difference),zeroM.N),Vr=fields.Vr,Vz=fields.Vz,tail=v.tail)
    end
    Tr=functional_bump(r,Ir,f.r_spec);Uz=functional_bump(z,Iz,f.z_spec)
    # Outside the vertical support all fields, including the scalar primitive, vanish.
    if sup(Iz)≤inf(f.z_spec[1]) || inf(Iz)≥sup(f.z_spec[4])
        return (q=zeroM,qnum=zeroM,δq=zeroM,Vr=zeroM,Vz=zeroM,tail=zeroM)
    end
    br=bases[f.r_basis];bz=bases[f.z_basis]
    radial_inside=!(sup(Ir)≤inf(f.r_spec[1]) || inf(Ir)≥sup(f.r_spec[4]))
    if f.kind=="vector" && !radial_inside
        return (q=zeroM,qnum=zeroM,δq=zeroM,Vr=zeroM,Vz=zeroM,tail=zeroM)
    end
    function mode(b,m,X,I)
        key=(b.name,objectid(b.modes),m)
        get!(cache,key) do
            functional_mode(b,m,X,I;N=N)
        end
    end
    function component(sigma,xi,which)
        value=zeroM;dr=zeroM;dz=zeroM;primitive=zeroM;tail=zeroM
        for n in eachindex(xi)
            m,k=sigma[n,1],sigma[n,2]
            (!radial_inside && m!=0) && continue
            Z=getproperty(mode(bz,k,z,Iz),which)
            c=interval(Float64,xi[n])
            if radial_inside
                R=getproperty(mode(br,m,r,Ir),which)
                value=value+c*(R.f*Z.f)
                dr=dr+c*(R.d*Z.f)
                dz=dz+c*(R.f*Z.d)
                if f.kind=="scalar"
                    F=get!(cache,(br.name,objectid(br.modes),m,which,:primitive)) do
                        m==0 ? br.constant*(functional_radial_primitive(r,Ir,f.r_spec)*r^I"-3") : -(I"1"/br.modes[m].lambda)*(Tr.f*R.d)
                    end
                    primitive=primitive+c*(F*Z.f)
                end
            elseif f.kind=="scalar" && inf(Ir)≥sup(f.r_spec[4]) && m==0
                H=br.constant*functional_bump_integral(r,f.r_spec,4)
                primitive=primitive+c*(H*r^I"-3"*Z.f)
            end
            if f.kind=="scalar" && m==0
                H=br.constant*functional_bump_integral(r,f.r_spec,4)
                tail=tail+c*(H*Z.f)
            end
        end
        return (f=value,d=dr,z=dz,primitive=primitive,tail=tail)
    end
    # Recompute the exact mathematical scaling, rather than treating a serialized
    # prefactor interval as a newly chosen constant. On a reflected positive-z cell
    # kappa=sqrt(2) and the parity projection supplies the factor1/2.
    ν=interval(Float64,f.nu);πb=RESPONSE_PI64
    factor=f.symmetric_cell ? sqrt(ν)/(sqrt(I"2")*πb) : sqrt(ν)/(I"2"*πb)
    function assemble(which)
        X1=component(f.sigma1,f.xi1,which)
        if f.kind=="vector"
            if !radial_inside
                return (q=zeroM,Vr=zeroM,Vz=zeroM,tail=zeroM)
            end
            X2=component(f.sigma2,f.xi2,which)
            V1=factor*(Tr.f*Uz.f*X1.f);V2=factor*(Tr.f*Uz.f*X2.f)
            q=-factor*(Uz.f*(Tr.d*X1.f+Tr.f*X1.d+I"3"*(Tr.f*X1.f)*r^I"-1")+Tr.f*(Uz.d*X2.f+Uz.f*X2.z))
            return (q=q,Vr=V1,Vz=V2,tail=zeroM)
        elseif f.kind=="scalar"
            return (q=factor*(Tr.f*Uz.f*X1.f),Vr=factor*(Uz.f*X1.primitive),Vz=zeroM,tail=factor*(Uz.f*X1.tail))
        end
        error("unknown functional kind $(f.kind)")
    end
    num=assemble(:num);exact=assemble(:exact)
    difference=exact.q-num.q
    δq=Model(zeroM.Π,sup_bound(difference),zeroM.N)
    return (q=exact.q,qnum=num.q,δq=δq,Vr=exact.Vr,Vz=exact.Vz,tail=exact.tail)
end
