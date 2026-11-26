if !isdefined(@__MODULE__,:PROOF_DATA)
    const PROOF_DATA=get(ENV,"PROOF_DATA",normpath(joinpath(@__DIR__,"..","data")))
    const COERCIVITY_RESULTS=get(ENV,"COERCIVITY_RESULTS",joinpath(@__DIR__,"results"))
    ENV["UNWEIGHTED_DATA"]=PROOF_DATA
    ENV["RESPONSE_DATA_FOLDER"]=PROOF_DATA
end
