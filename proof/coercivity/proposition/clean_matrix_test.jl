# Optional cross-check. Run with the environment of a clean run after Boxes.jld2
# has been assembled. This test computes one representative rectangle per row in
# memory at two independent Taylor contractions and never writes production data.
include("clean_matrices.jl")

function test_matrix_rectangles()
    inputs,boxes,_=clean_matrix_inputs()
    weights=clean_source_weights(inputs.coercivity.rows)
    length(weights)==58 || error("incomplete fixed weight table")
    for row in inputs.coercivity.rows
        ctx=clean_geometry(inputs,boxes,row)
        parent=cld(length(ctx.radial),2)
        u=ctx.radial[parent][1:2]
        v=ctx.axial[parent][1:2]
        cell=merge(ctx,(U=u,base_V=v,radial=[u],axial=[v]))
        moments=prepare_moments(inputs,cell,[1])
        G,G2=matrix_strip(inputs,boxes,cell,1,moments)
        baseline=gram_column(inputs.profile,row,inputs.coercivity.bases,boxes,inputs.far,u,v,1)
        for (fresh,reference) in ((G,baseline.G),(G2,baseline.G2))
            all(clean_finite,intersect_interval.(fresh,reference)) ||
                error("independent moment contractions disagree in $(row.id)")
            unpacked=matrix_unpack(matrix_pack(fresh),ctx.K)
            all(isequal_interval.(unpacked,fresh)) || error("symmetric packing changed the matrix")
        end
        println(row.id," independent rectangle check passed (",ctx.K," modes)")
    end
end
abspath(PROGRAM_FILE)==(@__FILE__) && test_matrix_rectangles()
