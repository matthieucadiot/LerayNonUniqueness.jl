# Run once to save outward Float64 interval coefficients for the coercivity workers.
include(joinpath(@__DIR__,"table2","Response_setup.jl"))
setprecision(128)

function check_conversion(A,B)
    size(A)==size(B) || error("the compact coefficient dimensions changed")
    all(i->isguaranteed(B[i])&&isfinite(inf(B[i]))&&isfinite(sup(B[i]))&&
        inf(B[i])<=inf(A[i])<=sup(A[i])<=sup(B[i]),eachindex(A)) || error("a Float64 interval lost its original coefficient")
    return length(A)
end

function check_profile_conversion(p,q)
    n=0
    for k in (:Ψ,:Ω,:Ψ̃,:Ω̃)
        hasproperty(p,k) || continue
        n+=check_conversion(getproperty(p,k),getproperty(q,k))
    end
    for k in keys(p)
        k in (:Ψ,:Ω,:Ψ̃,:Ω̃) || getproperty(p,k)===getproperty(q,k) || error("profile metadata changed")
    end
    check_compact_profile(q)
    return n
end

function save_compact(file,r,source,n)
    r=merge(Dict{String,Any}(r),Dict("schema"=>"coercivity-float64-v1","coefficient_precision"=>53,
        "topology_precision"=>128,"conversion_checked"=>true,"source_sha256"=>file_hash(source),
        "conversion_signature"=>file_hash(joinpath(@__DIR__,"Compact_data.jl")),"converted_coefficients"=>n))
    mkpath(dirname(file));temporary,io=mktemp(dirname(file));close(io)
    try
        jldopen(temporary,"w";compress=true) do output
            for (key,value) in r
                output[key]=value
            end
        end
        mv(temporary,file;force=true)
    finally
        isfile(temporary)&&rm(temporary)
    end
    println(basename(file),": ",n," outward coefficient enclosures, ",filesize(file)," bytes")
end

function main()
    source=joinpath(PROOF_DATA,"Profile","Interval_data.jld2")
    p=load_data(dirname(source));q=compact_profile(p)
    n=check_profile_conversion(p,q)
    save_compact(joinpath(COMPACT_DATA,"Profile_float64.jld2"),Dict("data"=>q),source,n)
    p=nothing;q=nothing;GC.gc()

    source=joinpath(PROOF_DATA,"Proposition_4_21_inputs.jld2")
    r=load(source);p=r["profile"];q=compact_profile(p)
    n=check_profile_conversion(p,q)
    bases=compact_bases(r["bases"])
    for (name,b) in r["bases"],(i,m) in enumerate(b.modes)
        target=bases[name].modes[i]
        for k in (:phi_a,:psi_a,:psi_b,:phi_b)
            n+=check_conversion(getproperty(m,k),getproperty(target,k))
        end
        for k in keys(m)
            k in (:phi_a,:psi_a,:psi_b,:phi_b) || getproperty(m,k)===getproperty(target,k) || error("spectral seed or metadata changed")
        end
    end
    r["profile"]=q;r["bases"]=bases
    save_compact(joinpath(COMPACT_DATA,"Proposition_float64.jld2"),r,source,n)
end
main()
