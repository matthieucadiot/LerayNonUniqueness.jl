# Pointwise coefficient checks shared by the column drivers.

function verify_farfield_cell(data,box,far,depth=0)
    F=coefficient_models(data,box,12,far)
    sup(functional_range(F.R2))<40000 && return true
    check_farfield_box(F) && return true
    depth<16 || error("far-field inequalities not established on $box")
    C,E=bisect_box(box)
    verify_farfield_cell(data,C,far,depth+1)
    verify_farfield_cell(data,E,far,depth+1)
    return true
end


# Cover the physical rectangle by OUTWARD transformed endpoints. In particular,
# keep the thin slabs between the two rounded endpoints of an irrational cut.
# The topology remains BigFloat; a rounded Float64 physical-range overlap is
# not a valid reason to impose a support inequality on an entire adjacent cell.
function bounded_domain_box(rectangle)
    r0,r1,z0,z1=interval.(BigFloat,rectangle)
    all(isguaranteed,(r0,r1,z0,z1)) || error("unguaranteed bounded-check domain")
    inf(r0)≥0 && inf(z0)≥0 && sup(r0)<inf(r1) && sup(z0)<inf(z1) ||
        error("invalid bounded-check domain")
    a,b,c,d=coordinate.((r0,r1,z0,z1))
    return interval.(BigFloat,(inf(a),sup(b),inf(c),sup(d)))
end

function bounded_domain_intersection(box,domain)
    all(x->isguaranteed(x)&&isthin(x),box) || error("bounded cell endpoints must be thin")
    all(x->isguaranteed(x)&&isthin(x),domain) || error("bounded domain endpoints must be thin")
    a=max(inf(box[1]),inf(domain[1]));b=min(sup(box[2]),sup(domain[2]))
    c=max(inf(box[3]),inf(domain[3]));d=min(sup(box[4]),sup(domain[4]))
    # Checks are made on CLOSED positive-area rectangles. A zero-area contact
    # is already covered by the closed neighboring intersection on the required
    # side of the boundary in the complete bounds grid. Do not construct a
    # zero-width Taylor cell, which would divide by hu=0 or hv=0.
    a<b && c<d || return nothing
    return interval.(BigFloat,(a,b,c,d))
end

function bounded_cell_tasks(inputs,box)
    D=inputs.coercivity
    tasks=Tuple{NamedTuple,Interval,NTuple{4,Interval{BigFloat}}}[]
    function add_task(row,w,rectangle)
        required=bounded_domain_intersection(box,bounded_domain_box(rectangle))
        required===nothing || push!(tasks,(row,w,required))
    end
    # These are the same exterior rectangles as verify_confinement. Their
    # inequalities intentionally apply outside the support rectangles below.
    odd=merge(first(D.rows),(η=I"0",))
    scalar=merge(first(filter(row->row.s==1,D.rows)),(η=I"0",))
    for rectangle in ((0,1,0,200),(42,200,0,200),(1,42,20,200))
        add_task(odd,I"0",rectangle)
    end
    for rectangle in ((0,I"0.2",0,200),(152,200,0,200),(I"0.2",152,35,200))
        add_task(scalar,I"0",rectangle)
    end
    for row in D.rows
        r=D.bases[row.r_basis].nodes;z=D.bases[row.z_basis].nodes
        add_task(row,row.W,(r[1],r[end],max(big(I"0"),z[1]),z[end]))
    end
    return tasks
end

function verify_bounded_cell(inputs,B,box)
    p=inputs.profile;far=inputs.far
    tasks=bounded_cell_tasks(inputs,box)
    isempty(tasks) && return true
    first_box=first(tasks)[3]
    fields=Dict(first_box=>coefficient_models(p,first_box,12,far))
    zetas=Dict{typeof(first_box),Interval{Float64}}()
    for (row,w,required) in tasks
        F=get!(fields,required) do
            coefficient_models(p,required,12,far)
        end
        ζ=row.block==:odd ? get!(() -> zeta_bound(B,required),zetas,required) : I"0"
        M=matrix_field(F,row,ζ;positivity=true)
        positive_matrix(M,w) && continue
        verify_positive_box(p,required,far,row,B,w)
    end
    return true
end

