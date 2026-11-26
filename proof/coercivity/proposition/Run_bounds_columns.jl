# Shared coefficient evaluation for (iii) and all fifteen bounds (iv).
# julia --project=. proposition/Run_bounds_columns.jl FIRST LAST
isdefined(@__MODULE__,:load_proof_inputs) || include("Setup.jl")
isdefined(@__MODULE__,:read_proof) || include("Files.jl")


deadline=proof_deadline()
inputs=load_proof_inputs();B=read_proof(joinpath(OUTPUT,"Boxes.jld2"))["boxes"]
U,V=bounds_grid(inputs,B)
first_column=isempty(ARGS) ? 1 : parse(Int,ARGS[1])
last_column=length(ARGS)<2 ? length(U)-1 : parse(Int,ARGS[2])
1≤first_column≤last_column<length(U) || error("invalid bounded-column range")
Threads.@threads for i in first_column:last_column
    file=joinpath(OUTPUT,"bounds_columns","column_$(i).jld2")
    if isfile(file)
        r=read_proof(file)
        r["column"]==i && r["rows"]==collect(1:length(V)-1) && r["verified"] || error("incomplete bounded column")
    else
        time()<deadline || continue
        for j in 1:length(V)-1
            verify_bounded_cell(inputs,B,(U[i],U[i+1],V[j],V[j+1]))
        end
        save_proof(file;column=i,rows=collect(1:length(V)-1),verified=true)
    end
    println("(iii),(iv) bounded column ",i," / ",length(U)-1);flush(stdout)
end

all(i->isfile(joinpath(OUTPUT,"bounds_columns","column_$(i).jld2")),first_column:last_column) || exit(75)
