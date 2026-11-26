# Exact subdivisions of the original response cells at spectral and support junctions.
# Every refinement is an exact partition of the existing transformed-coordinate cell.

function refinement_validate_quantities(S; odd=true)
    length(S)==8 || error("expected all eight Table 2 quantities")
    all(x -> x isa Interval{Float64} && isguaranteed(x) &&
        isfinite(inf(x)) && isfinite(sup(x)) && inf(x)<=sup(x), S) ||
        error("invalid or nonguaranteed Float64 cell quantity")
    all(k -> (!odd && k==6) || inf(S[k])>=0,eachindex(S)) ||
        error("a squared norm is negative")
    return S
end

function refinement_check_partition(original,parts)
    length(original)==2 && length(parts)>=2 || error("invalid partition dimensions")
    all(x -> x isa Interval{BigFloat} && isguaranteed(x) && isthin(x) &&
        isfinite(inf(x)), parts) || error("partition endpoints must be finite guaranteed thin BigFloat intervals")
    isequal(parts[1],original[1]) && isequal(parts[end],original[end]) ||
        error("refinement changed an original outer endpoint")
    all(i -> sup(parts[i])<inf(parts[i+1]),1:length(parts)-1) ||
        error("partition is not strictly increasing")
    # Adjacent children are represented by consecutive entries of this ONE vector.
    # Thus they share the identical exact point, with no gaps or missing strip.
    return true
end

function refinement_junction(node)
    x=abs(interval(BigFloat,node))
    isguaranteed(x) && isfinite(inf(x)) && isfinite(sup(x)) || error("invalid physical junction")
    return (interval(BigFloat,1)+interval(BigFloat,12)*x*x)^interval(BigFloat,1//4)-interval(BigFloat,1)
end

function refinement_coordinate_partition(original,junctions; ulps=16)
    ulps isa Integer && ulps>0 || error("collar width must be a positive integer")
    refinement_check_partition(original,original)
    a,b=inf(original[1]),inf(original[2])
    cuts=BigFloat[]; touched=0; unresolved=0
    for J in junctions
        J isa Interval{BigFloat} && isguaranteed(J) && inf(J)<=sup(J) || error("invalid transformed junction")
        # No subnormal collar is inserted at the origin. The original axis cell
        # and the production axis chart are retained intact there.
        sup(J)<=0 && continue
        lo=prevfloat(Float64(inf(J)),ulps)
        hi=nextfloat(Float64(sup(J)),ulps)
        isfinite(lo) && isfinite(hi) || error("nonfinite collar")
        left=BigFloat(lo);right=BigFloat(hi)
        max(a,left)<=min(b,right) || continue
        touched+=1
        # Include the outward junction enclosure as well as its two collars.
        # If the original cell is microscopic, there may be no interior cut;
        # retaining it whole still covers the domain rigorously.
        n=length(cuts)
        for y in (left,inf(J),sup(J),right)
            a<y<b && push!(cuts,y)
        end
        length(cuts)==n && (unresolved+=1)
    end
    points=sort(unique(cuts))
    parts=vcat([original[1]],interval.(BigFloat,points),[original[2]])
    refinement_check_partition(original,parts)
    return (parts=parts,touched=touched,unresolved=unresolved)
end

function refinement_axis_plan(grid,physical_nodes; ulps=16)
    length(grid)>=2 || error("empty grid")
    junctions=[refinement_junction(x) for x in physical_nodes]
    # Identical reflected nodes need not create duplicate reports/cuts.
    junctions=unique(x -> (inf(x),sup(x)),junctions)
    partitions=Vector{Vector{Interval{BigFloat}}}()
    touched=Int[];refined=Int[];microscopic=Int[]
    for i in 1:length(grid)-1
        result=refinement_coordinate_partition(grid[i:i+1],junctions;ulps=ulps)
        push!(partitions,result.parts)
        result.touched>0 && push!(touched,i)
        length(result.parts)>2 && push!(refined,i)
        result.touched>0 && length(result.parts)==2 && push!(microscopic,i)
    end
    records=[Dict("physical_lower"=>string(inf(abs(interval(BigFloat,x)))),
                  "physical_upper"=>string(sup(abs(interval(BigFloat,x)))),
                  "transformed_lower"=>string(inf(refinement_junction(x))),
                  "transformed_upper"=>string(sup(refinement_junction(x)))) for x in physical_nodes]
    return (partitions=partitions,touched=touched,refined=refined,microscopic=microscopic,junctions=records)
end

function refinement_nodes(f,bases,axis)
    axis in (:u,:v) || error("unknown coordinate axis")
    spec=axis==:u ? f.r_spec : f.z_spec
    basis=bases[axis==:u ? f.r_basis : f.z_basis]
    lower=max(BigFloat(0),BigFloat(inf(spec[1])))
    upper=BigFloat(sup(spec[4]))
    result=Interval{BigFloat}[]
    for value in vcat(collect(basis.nodes),collect(spec))
        x=abs(interval(BigFloat,value))
        sup(x)>0 && max(inf(x),lower)<=min(sup(x),upper) && push!(result,x)
    end
    return unique(x -> (inf(x),sup(x)),result)
end

function refinement_physical_interval(endpoints)
    return interval(inf(sqrt(max(interval(BigFloat,0),U0(endpoints[1])))),sup(sqrt(U0(endpoints[end]))))
end

function refinement_prepare_case(t,par,case,plan)
    odd=case.odd;index=case.index;col=case.column
    data=response_profile(t,odd);rd=odd ? t.odd : t.even;qs=odd ? t.qo : t.qe
    1<=col<length(data.U) && 1<=index<=rd.n || error("invalid validation case")
    uparts=plan.U.partitions[col]
    d=merge(data,(U=uparts,))
    C=[column_data(k,d,par) for k in 1:length(uparts)-1]
    CR=[response_basis_column(k,rd,uparts) for k in 1:length(uparts)-1]
    far=far_field_data(par.prec,data.λ̃)
    return (case=case,par=par,data=data,rd=rd,qs=qs,f=qs.functionals[index],far=far,
        ext=response_farfield(index,rd,far),uparts=uparts,C=C,CR=CR,
        vparts=plan.V.partitions,nv=length(data.V)-1)
end

function refinement_evaluate_row(ctx,j)
    1<=j<=ctx.nv || error("invalid original row")
    started=time();vp=ctx.vparts[j]
    refinement_check_partition(ctx.data.U[ctx.case.column:ctx.case.column+1],ctx.uparts)
    refinement_check_partition(ctx.data.V[j:j+1],vp)
    d=merge(ctx.data,(U=ctx.uparts,V=vp))
    total=zeros(Interval{Float64},8);children=0
    for a in eachindex(ctx.C),b in 1:length(vp)-1
        C=ctx.C[a];CR=ctx.CR[a]
        R=row_data(b,d,ctx.par);RR=response_basis_row(b,ctx.rd,vp)
        P=profile_fields_on_cell(d,ctx.far,C,R,ctx.par)
        Q=functional_on_cell(ctx.f,ctx.qs.eigenfunctions,C,R,ctx.par)
        S=vcat(response_cell_quantities(ctx.case.index,ctx.rd,CR,RR,C,R,P,ctx.ext,Q.q,ctx.par),
               functional_cell_integrals(Q,ctx.case.odd,C,R))
        refinement_validate_quantities(S;odd=ctx.case.odd)
        total=total+S;children+=1
    end
    children==(length(ctx.uparts)-1)*(length(vp)-1) || error("a child cell was omitted")
    refinement_validate_quantities(total;odd=ctx.case.odd)
    return (S=total,children=children,coverage_checked=true,elapsed_seconds=time()-started)
end

function refinement_flatten_partition(grid,parts)
    length(parts)==length(grid)-1 || error("partition grid size mismatch")
    result=[grid[1]]
    for j in eachindex(parts)
        refinement_check_partition(grid[j:j+1],parts[j])
        isequal(result[end],parts[j][1]) || error("gap between original cells")
        append!(result,parts[j][2:end])
    end
    refinement_check_partition([grid[1],grid[end]],result)
    return result
end

function refinement_tail_rows(t,par,index,plan)
    f=t.qe.functionals[index];f.kind=="scalar" || error("only scalar tails requested")
    data=response_profile(t,false)
    refined_grid=refinement_flatten_partition(data.V,plan.V.partitions)
    d=merge(data,(V=refined_grid,))
    original=[row_data(j,data,par) for j in 1:length(data.V)-1]
    refined=[row_data(j,d,par) for j in 1:length(refined_grid)-1]
    return (original=original,refined=refined)
end
