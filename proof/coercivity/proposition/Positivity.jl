# Proposition 4.21(iii),(iv): interval diagonal/determinant tests.

function positive_matrix(M,w=I"0")
    if length(M)==1
        q=functional_range(M[1]+w)
        isguaranteed(q) || error("unguaranteed scalar range")
        return inf(q)≥0
    end
    a,b,c=M; a=a+w; c=c+w
    ranges=(functional_range(a),functional_range(c),functional_range(a*c-b*b))
    all(isguaranteed,ranges) || error("unguaranteed matrix range")
    return all(q->inf(q)≥0,ranges)
end

function verify_positive_box(data,box,far,row,B,w=I"0",depth=0)
    F=coefficient_models(data,box,12,far)
    ζ=row.block==:odd ? zeta_bound(B,box) : I"0"
    M=matrix_field(F,row,ζ;positivity=true)
    positive_matrix(M,w) && return true
    depth<16 || error("positivity not established for $(row.id) on $box")
    C,D=bisect_box(box)
    verify_positive_box(data,C,far,row,B,w,depth+1)
    verify_positive_box(data,D,far,row,B,w,depth+1)
    return true
end

# Refine at spline knots, box jumps, and outward endpoints of physical boundaries.
function proof_grid(grid,lo,hi,extra=())
    a=coordinate(interval(BigFloat,lo)); b=coordinate(interval(BigFloat,hi))
    x=BigFloat[inf(a),sup(a),inf(b),sup(b)]
    append!(x,inf.(filter(t->inf(a)<inf(t)<sup(b),grid)))
    for v in extra
        q=coordinate(interval(BigFloat,v))
        for e in (inf(q),sup(q))
            inf(a)<e<sup(b) && push!(x,e)
        end
    end
    return interval.(BigFloat,sort(unique(x)))
end

function verify_rectangle(data,far,row,B,r0,r1,z0,z1,w)
    U=proof_grid(data.U,r0,r1)
    V=proof_grid(data.V,z0,z1)
    # Every zeta jump is already a knot of B.V; B.U adds the radial box partition.
    append!(U,filter(t->inf(first(U))<inf(t)<sup(last(U)),B.U));sort!(U,by=inf);unique!(U)
    append!(V,filter(t->inf(first(V))<inf(t)<sup(last(V)),B.V));sort!(V,by=inf);unique!(V)
    for i in 1:length(U)-1,j in 1:length(V)-1
        verify_positive_box(data,(U[i],U[i+1],V[j],V[j+1]),far,row,B,w)
    end
    return true
end

function verify_confinement(data,far,rows,B)
    odd=merge(first(rows),(η=I"0",))
    scalar=merge(first(filter(r->r.s==1,rows)),(η=I"0",))
    for (row,rectangles) in ((odd,((0,1,0,200),(42,200,0,200),(1,42,20,200))),
                            (scalar,((0,I"0.2",0,200),(152,200,0,200),(I"0.2",152,35,200))))
        for rectangle in rectangles
            verify_rectangle(data,far,row,B,rectangle...,I"0")
        end
    end
    return true
end

function verify_support_bounds(data,far,rows,bases,B)
    for row in rows
        r=bases[row.r_basis].nodes; z=bases[row.z_basis].nodes
        z0=max(I"0",z[1])
        verify_rectangle(data,far,row,B,r[1],r[end],z0,z[end],row.W)
        println("(iv) ",row.id," verified")
    end
    return true
end
