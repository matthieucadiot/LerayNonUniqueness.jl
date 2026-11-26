# julia --project=. proposition/Run_boxes.jl [first_column last_column]
isdefined(@__MODULE__,:load_proof_inputs) || include("Setup.jl")
isdefined(@__MODULE__,:read_proof) || include("Files.jl")

deadline=proof_deadline()
inputs=load_proof_inputs()
data=inputs.profile;far=inputs.far
U=box_partition(data.U,2482); V=box_partition(data.V)
first_column=isempty(ARGS) ? 1 : parse(Int,ARGS[1])
last_column=length(ARGS)<2 ? length(U)-1 : parse(Int,ARGS[2])
1≤first_column≤last_column≤length(U)-1 || error("invalid column range")
folder=joinpath(OUTPUT,"box_columns");mkpath(folder)
Threads.@threads for i in first_column:last_column
    file=joinpath(folder,"boxes_$(i).jld2")
    if isfile(file)
        d=read_proof(file)
        d["i"]==i && all(isequal_interval.(d["U"],U)) && all(isequal_interval.(d["V"],V)) || error("wrong saved box column")
        verify_box_column(d["lower"],d["box"])
    else
        time()<deadline || continue
        C=compute_box_column(data,U,V,i,far)
        save_proof(file;i=i,U=U,V=V,lower=C.lower,box=C.box)
    end
    println("box column ",i," / ",length(U)-1)
    flush(stdout)
end
all(i->isfile(joinpath(folder,"boxes_$(i).jld2")),first_column:last_column) || exit(75)
