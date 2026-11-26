# Proposition 4.21(ii). Box boundaries are defined by thin u,v endpoints;
# their physical values are r=sqrt(U₀(u)), z=sqrt(U₀(v)).

function box_partition(grid,K=0)
    b=interval(BigFloat,sup(coordinate(I"200")))
    x=sort(unique(vcat(inf.(filter(t->sup(t)<inf(b),grid)),[inf(b)])))
    while K>0 && length(x)<K+1
        j=argmax(diff(x)); insert!(x,j+1,(x[j]+x[j+1])/2)
    end
    K==0 || length(x)==K+1 || error("more than $K original knot intervals below r=200")
    return interval.(BigFloat,x)
end

# A midpoint is only a proposed new boundary, never a bound on a function.
function bisect_box(box)
    a,b,c,d=box
    if sup(b-a)≥sup(d-c)
        m=interval(BigFloat,mid((a+b)/I"2"))
        return ((a,m,c,d),(m,b,c,d))
    end
    m=interval(BigFloat,mid((c+d)/I"2"))
    return ((a,b,c,m),(a,b,m,d))
end

function e0_lower_bound(data,box,far,N=12,depth=0)
    F=coefficient_models(data,box,N,far)
    q=functional_range(F.e0)
    isguaranteed(q) && isfinite(inf(q)) || error("invalid enclosure of e₀")
    # In the positive exterior, subdivision prevents a spurious negative box reaching 200.
    if inf(q)<0 && sup(q)>0 && inf(U0(box[3]))>1 && depth<4
        B,C=bisect_box(box)
        return min(e0_lower_bound(data,B,far,N,depth+1),e0_lower_bound(data,C,far,N,depth+1))
    end
    return interval(Float64,inf(q))
end

# Two heights: h₂ bounds the outer part, h₁+h₂ the inner part.
# Minimize h₁L₁²+h₂L₂² over the available vertical grid levels.
function two_boxes(lower,V)
    n=[max(I"0",-x) for x in lower]
    last=findlast(x->sup(x)>0,n)
    last===nothing && return (j1=0,j2=0,h1=I"0",h2=I"0",L1=I"0",L2=I"0")
    last<length(n) || error("the enclosure of the negative part reaches z=200; refine the cells")
    L2=interval(Float64,sqrt(U0(V[last+1]))); sup(L2)<200 || error("L₂ must be below 200")
    best=Float64(Inf); result=nothing
    for j in 0:last
        h2=j==last ? I"0" : interval(Float64,nextfloat(maximum(sup.(n[j+1:last]))))
        h1=j==0 ? I"0" : interval(Float64,nextfloat(max(Float64(0),sup(maximum(n[1:j])-h2))))
        L1=j==0 ? I"0" : interval(Float64,sqrt(U0(V[j+1])))
        cost=sup(h1*L1^2+h2*L2^2)
        if cost<best
            best=cost; result=(j1=j,j2=last,h1=h1,h2=h2,L1=L1,L2=L2)
        end
    end
    verify_box_column(lower,result)
    return result
end

function verify_box_column(lower,b)
    0≤b.j1≤b.j2≤length(lower) || error("invalid box levels")
    all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),lower) || error("invalid e₀ lower bounds")
    all(x->isguaranteed(x)&&isfinite(sup(x))&&inf(x)≥0,(b.h1,b.h2,b.L1,b.L2)) || error("invalid box constants")
    sup(b.L1)≤sup(b.L2)<200 || error("invalid physical box levels")
    for j in eachindex(lower)
        h=j≤b.j1 ? b.h1+b.h2 : j≤b.j2 ? b.h2 : I"0"
        inf(lower[j]+h)≥0 || error("box inequality failed in vertical cell $j")
    end
    return true
end

function compute_box_column(data,U,V,i,far)
    lower=[e0_lower_bound(data,(U[i],U[i+1],V[j],V[j+1]),far) for j in 1:length(V)-1]
    return (lower=lower,box=two_boxes(lower,V))
end

function assemble_boxes(folder,U,V)
    columns=map(1:length(U)-1) do i
        c=read_proof(joinpath(folder,"boxes_$(i).jld2"))
        c["i"]==i && isequal_interval.(c["U"],U)==fill(true,length(U)) &&
            isequal_interval.(c["V"],V)==fill(true,length(V)) || error("box column $i has a different grid")
        q=c["box"]
        length(c["lower"])==length(V)-1 || error("box column $i has incomplete vertical coverage")
        verify_box_column(c["lower"],q)
        for (j,L) in ((q.j1,q.L1),(q.j2,q.L2))
            expected=j==0 ? I"0" : interval(Float64,sqrt(U0(V[j+1])))
            isequal_interval(L,expected) || error("a saved box level has changed")
        end
        return q
    end
    r=vcat(sqrt.(U0.(U[1:end-1])),[I"200"])
    B=(r=r,U=U,V=V,columns=columns)
    length(B.columns)==2482 || error("Proposition 4.21 specifies 2482 boxes")
    return B
end

# An interval upper bound also works for cells crossing a jump. Gram cells are
# refined at all jumps, where this routine returns the defining constant.
function zeta_bound(B,box)
    a,b,c,d=box; value=I"0"
    for i in eachindex(B.columns)
        inf(b)>inf(B.U[i]) && sup(a)<sup(B.U[i+1]) || continue
        q=B.columns[i]; h=I"0"
        q.j1>0 && inf(c)<sup(B.V[q.j1+1]) && (h=h+proof_float64(q.h1)*proof_float64(q.L1)^2)
        q.j2>0 && inf(c)<sup(B.V[q.j2+1]) && (h=h+proof_float64(q.h2)*proof_float64(q.L2)^2)
        value=max(value,h)
    end
    return I"4"*value/interval(Float64,π)^2
end
