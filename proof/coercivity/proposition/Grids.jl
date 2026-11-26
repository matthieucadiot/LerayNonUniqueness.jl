# Thin transformed endpoints, refined at every spline knot and physical cut.
function matrix_grid(grid,lo,hi,nodes,boxgrid)
    U=proof_grid(grid,lo,hi,nodes)
    append!(U,filter(t->inf(first(U))<inf(t)<sup(last(U)),boxgrid))
    return interval.(BigFloat,sort(unique(inf.(U))))
end

function bounds_grid(inputs,B)
    r=interval.(BigFloat,[0,1,42,152,200]);push!(r,big(I"0.2"))
    z=interval.(BigFloat,[0,20,35,200])
    for b in values(inputs.coercivity.bases)
        append!(b.radial ? r : z,abs.(b.nodes))
    end
    U=matrix_grid(inputs.profile.U,big(I"0"),big(I"200"),r,B.U)
    V=matrix_grid(inputs.profile.V,big(I"0"),big(I"200"),z,B.V)
    return U,V
end

function matrix_grids(inputs,B,row)
    br=inputs.coercivity.bases[row.r_basis];bz=inputs.coercivity.bases[row.z_basis]
    U=matrix_grid(inputs.profile.U,br.nodes[1],br.nodes[end],br.nodes,B.U)
    V=matrix_grid(inputs.profile.V,max(big(I"0"),bz.nodes[1]),bz.nodes[end],bz.nodes,B.V)
    return U,V
end

function check_column_coverage(columns,n)
    counts=zeros(Int,n)
    for group in columns, i in group
        i isa Integer && 1≤i≤n || error("invalid column index")
        counts[i]+=1
    end
    all(==(1),counts) || error("missing or duplicate columns")
    return true
end
