# julia --project=. proposition/Run_farfield_columns.jl FIRST LAST
isdefined(@__MODULE__,:load_proof_inputs) || include("Setup.jl")
isdefined(@__MODULE__,:read_proof) || include("Files.jl")


deadline=proof_deadline()
inputs=load_proof_inputs();p=inputs.profile
first_column=isempty(ARGS) ? 1 : parse(Int,ARGS[1])
last_column=length(ARGS)<2 ? length(p.U)-1 : parse(Int,ARGS[2])
1≤first_column≤last_column<length(p.U) || error("invalid far-field column range")
Threads.@threads for i in first_column:last_column
    file=joinpath(OUTPUT,"farfield_columns","column_$(i).jld2")
    if isfile(file)
        r=read_proof(file)
        r["column"]==i && r["rows"]==collect(1:length(p.V)-1) && r["verified"] || error("incomplete far-field column")
    else
        time()<deadline || continue
        for j in 1:length(p.V)-1
            sup(U0(p.U[i+1])+U0(p.V[j+1]))≥40000 || continue
            verify_farfield_cell(p,(p.U[i],p.U[i+1],p.V[j],p.V[j+1]),inputs.far)
        end
        save_proof(file;column=i,rows=collect(1:length(p.V)-1),verified=true)
    end
    println("(i) far-field column ",i," / ",length(p.U)-1);flush(stdout)
end

all(i->isfile(joinpath(OUTPUT,"farfield_columns","column_$(i).jld2")),first_column:last_column) || exit(75)
