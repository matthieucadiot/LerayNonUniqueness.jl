# Compact enclosures for coercivity. The original residual and spectral data are retained.
using IntervalArithmetic,RadiiPolynomial,JLD2,SHA

const COMPACT_DATA=get(ENV,"COERCIVITY_DATA",normpath(joinpath(@__DIR__,"..","data","Coercivity")))
file_hash(file)=open(sha256,file) |> bytes2hex

function compact_profile(p)
    fields=(:Ψ,:Ω,:Ψ̃,:Ω̃)
    return merge(p,(; (k=>interval.(Float64,getproperty(p,k)) for k in fields if hasproperty(p,k))...))
end

function compact_bases(bases)
    return Dict(name=>merge(b,(modes=[merge(m,(;
        (k=>interval.(Float64,getproperty(m,k)) for k in (:phi_a,:psi_a,:psi_b,:phi_b))...))
        for m in b.modes],)) for (name,b) in bases)
end

function check_compact_profile(p)
    for k in (:Ψ,:Ω,:Ψ̃,:Ω̃)
        hasproperty(p,k) || continue
        A=getproperty(p,k)
        eltype(A)==Interval{Float64} && all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),A) ||
            error("invalid compact profile coefficients: $k")
        size(A)==(length(p.tu)-p.p-1,length(p.tv)-p.p-1) || error("invalid compact spline dimensions")
    end
    for k in (:tu,:tv,:U,:V)
        A=getproperty(p,k)
        all(x->x isa Interval{BigFloat}&&isguaranteed(x)&&isthin(x)&&precision(inf(x))==128,A) &&
            issorted(inf.(A)) || error("the original 128-bit topology was not retained")
    end
    return p
end

function check_compact_data(r,source)
    r["schema"]=="coercivity-float64-v1" && r["coefficient_precision"]==53 &&
        r["topology_precision"]==128 && r["conversion_checked"]===true || error("invalid compact data format")
    r["source_sha256"]==file_hash(source) || error("the compact data use different original inputs")
    r["conversion_signature"]==file_hash(@__FILE__) || error("rebuild the compact data after changing the conversion")
    return r
end

function load_coercivity_profile(folder)
    source=joinpath(folder,"Interval_data.jld2")
    file=joinpath(COMPACT_DATA,"Profile_float64.jld2")
    p=isfile(file) ? check_compact_data(load(file),source)["data"] : compact_profile(load_data(folder))
    return check_compact_profile(p)
end

function load_coercivity_proposition(source)
    file=joinpath(COMPACT_DATA,"Proposition_float64.jld2")
    if isfile(file)
        r=check_compact_data(load(file),source)
    else
        r=load(source)
        r["profile"]=compact_profile(r["profile"])
        r["bases"]=compact_bases(r["bases"])
    end
    check_compact_profile(r["profile"])
    for b in values(r["bases"]),m in b.modes,k in (:phi_a,:psi_a,:psi_b,:phi_b)
        A=getproperty(m,k)
        eltype(A)==Interval{Float64} && length(A)==201 &&
            all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),A) || error("invalid compact spectral coefficients")
    end
    return r
end
