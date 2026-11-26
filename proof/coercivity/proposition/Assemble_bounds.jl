isdefined(@__MODULE__,:load_proof_inputs) || include("Setup.jl")
isdefined(@__MODULE__,:read_proof) || include("Files.jl")
inputs=load_proof_inputs();p=inputs.profile
B=read_proof(joinpath(OUTPUT,"Boxes.jld2"))["boxes"]
U,V=bounds_grid(inputs,B)
for (kind,n,nv) in (("farfield",length(p.U)-1,length(p.V)-1),("bounds",length(U)-1,length(V)-1))
    columns=Vector{Int}[]
    for i in 1:n
        r=read_proof(joinpath(OUTPUT,kind*"_columns","column_$(i).jld2"))
        r["column"]==i && r["rows"]==collect(1:nv) && r["verified"] || error("incomplete $kind column $i")
        push!(columns,[r["column"]])
    end
    check_column_coverage(columns,n)
end
exterior=exterior_farfield(sqrt(U0(interval(BigFloat,p.L))))
exterior.ok || error("the analytic exterior estimate failed")
save_proof(joinpath(OUTPUT,"Bounds.jld2");farfield=exterior,finite_farfield=true,confinement=true,support_bounds=true)
println("(i),(iii),(iv) all columns verified")
