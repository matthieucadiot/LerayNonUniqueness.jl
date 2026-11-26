# Associate saved bounds with the formulas and input bytes that produced them.
function compute_proof_signature()
    files=filter(f->endswith(f,".jl") || endswith(f,".toml"),readdir(@__DIR__;join=true))
    append!(files,[joinpath(CORE,f) for f in readdir(CORE) if endswith(f,".jl")])
    append!(files,[joinpath(@__DIR__,"..",f) for f in ("Project.toml","Manifest.toml","Paths.jl")])
    push!(files,joinpath(PROOF_DATA,"Proposition_4_21_inputs.jld2"))
    push!(files,joinpath(@__DIR__,"..","Compact_data.jl"))
    cache=joinpath(COMPACT_DATA,"Proposition_float64.jld2")
    isfile(cache) && push!(files,cache)
    return Dict(relpath(f,dirname(@__DIR__))=>bytes2hex(sha256(read(f))) for f in files)
end

const SIGNATURE=compute_proof_signature()
proof_signature()=SIGNATURE

const PROOF_SAVE_LOCK=ReentrantLock()

function save_proof(file;kwargs...)
    lock(PROOF_SAVE_LOCK) do
        save_proof_unlocked(file;kwargs...)
    end
end

function save_proof_unlocked(file;kwargs...)
    mkpath(dirname(file))
    temporary,io=mktemp(dirname(file));close(io)
    try
        jldsave(temporary,true;arithmetic="Float64",input_prec=128,signature=proof_signature(),
                run_id=ENV["CLEAN_RUN_ID"],kwargs...)
        Base.Filesystem.rename(temporary,file)
    finally
        isfile(temporary) && rm(temporary)
    end
end

function read_proof(file)
    bytes=read(file)
    d=jldopen(bytes,"r") do f
        Dict{String,Any}(name=>f[name] for name in keys(f))
    end
    get(d,"arithmetic",nothing)=="Float64" && get(d,"input_prec",nothing)==128 || error("checkpoint arithmetic differs")
    get(d,"run_id",nothing)==ENV["CLEAN_RUN_ID"] || error("Proposition checkpoint belongs to a different run: $file")
    get(d,"signature",nothing)==proof_signature() || error("Proposition source or inputs differ: $file")
    return d
end

# Stop between complete columns, leaving the atomic checkpoints ready to resume.
function proof_deadline()
    seconds=parse(Float64,get(ENV,"PROOF_JOB_SECONDS","Inf"))
    seconds≥0 || error("invalid proposition soft time budget")
    return time()+seconds
end
