isdefined(@__MODULE__,:load_proof_inputs) || include("Setup.jl")
isdefined(@__MODULE__,:read_proof) || include("Files.jl")
p=load_proof_inputs().profile
U=box_partition(p.U,2482);V=box_partition(p.V)
folder=joinpath(OUTPUT,"box_columns")
boxes=assemble_boxes(folder,U,V)
save_proof(joinpath(OUTPUT,"Boxes.jld2");boxes=boxes)
println("(ii) All 2482 box columns verified")
