# Refine the compactified support at basis-piece junctions and enclose the monotone axial coordinate.
function strip_checked_points(original, points)
    all(x -> x isa Interval{BigFloat} && isguaranteed(x) && isthin(x) &&
             isfinite(inf(x)), points) || error("invalid integration cut")
    isequal_interval(first(original), first(points)) &&
        isequal_interval(last(original), last(points)) ||
        error("strip partition changed its outer endpoints")
    all(i -> sup(points[i]) < inf(points[i+1]), 1:length(points)-1) ||
        error("strip partition is not strictly increasing")
    return points
end

# Identical transformed coordinate and 16-binary64-ULP junction collars to
# the reviewed Table 2 repair. Every original parent is also bisected away from the exact axis when
# it has an interior representable midpoint. Keep every original grid point.
function strip_axis_plan(grid, physical_nodes)
    precision(BigFloat) == 128 || error("128-bit geometry required")
    strip_checked_points(grid, grid)
    junctions = [((interval(BigFloat,1) + interval(BigFloat,12) *
                   abs(interval(BigFloat,x))^2)^interval(BigFloat,1//4) -
                  interval(BigFloat,1)) for x in physical_nodes]
    partitions = Vector{Vector{eltype(grid)}}()
    near_junction = Bool[]
    for j in 1:length(grid)-1
        parent = grid[j:j+1]
        a,b = inf.(parent)
        midpoint = a + (b-a)/2
        cuts = BigFloat[a,b]
        !iszero(a) && a < midpoint < b && push!(cuts, midpoint)
        touched = false
        for z in junctions
            isguaranteed(z) && isfinite(inf(z)) && isfinite(sup(z)) ||
                error("invalid transformed junction")
            sup(z) <= 0 && continue  # Keep the axis chart at its exact endpoint.
            lower = BigFloat(prevfloat(Float64(inf(z)),16))
            upper = BigFloat(nextfloat(Float64(sup(z)),16))
            max(a,lower) <= min(b,upper) || continue
            touched = true
            for point in (lower, inf(z), sup(z), upper)
                a < point < b && push!(cuts,point)
            end
        end
        parts = interval.(BigFloat,sort!(unique!(cuts)))
        push!(partitions,strip_checked_points(parent,parts))
        push!(near_junction,touched)
    end
    refined = [grid[1]]
    for parts in partitions
        isequal_interval(last(refined), first(parts)) || error("gap between parents")
        append!(refined,parts[2:end])
    end
    strip_checked_points(grid,refined)
    return (; partitions,refined,near_junction)
end

function strip_validate_matrix(A, K)
    size(A) == (K,K) || error("unexpected strip Gram dimensions")
    all(x -> x isa Interval{Float64} && isguaranteed(x) && !isempty_interval(x) &&
             isfinite(inf(x)) && isfinite(sup(x)), A) ||
        error("invalid strip Gram entries")
    all(isequal_interval.(A,transpose(A))) || error("asymmetric strip Gram entries")
    return A
end

function strip_z_range(model,c,d)
    inf(c) >= 0 && sup(c) < inf(d) || error("invalid vertical slab")
    nonnegative = interval(BigFloat,0,Inf)
    lower = sqrt(intersect_interval(U0(c),nonnegative))
    upper = sqrt(intersect_interval(U0(d),nonnegative))
    physical = proof_float64(interval(BigFloat,inf(lower),sup(upper)))
    valid = intersect_interval(functional_range(model),physical)
    isguaranteed(valid) && !isempty_interval(valid) || error("inconsistent axial coordinates")
    return valid
end

function row_axial_plan(V, nodes)
    all_cuts = strip_axis_plan(V,nodes)
    partitions = [all_cuts.near_junction[j] ? all_cuts.partitions[j] : V[j:j+1]
                  for j in eachindex(all_cuts.partitions)]
    refined = [first(V)]
    for (j,parts) in pairs(partitions)
        strip_checked_points(V[j:j+1],parts)
        isequal_interval(last(refined),first(parts)) || error("gap before axial parent $j")
        append!(refined,parts[2:end])
    end
    strip_checked_points(V,refined)
    return (;partitions,refined,junctions=findall(all_cuts.near_junction))
end

