# Each odd coefficient is exactly (hi + lo)/denominator. Store both binary
# numerator limbs so close coefficients can cancel before the outward division.
using JLD2, SHA, IntervalArithmetic
const CLEAN_WEIGHTS_SHA = "5b3f70f66a57ce1ad934fe1b9bdf5a8a0130d4dbe19f9cb6b2034af4d309801b"
const CLEAN_INPUT_NAME = "Fixed_weighted_responses.jld2"
const CLEAN_BASE_FILES = ("functionals_odd.jld2","functionals_even.jld2","responses_odd.jld2","responses_even.jld2")
clean_hash(p) = bytes2hex(open(sha256,p))

function clean_design(folder)
    f = joinpath(folder,"Coercivity","selected_factors.tsv")
    clean_hash(f)==CLEAN_WEIGHTS_SHA || error("The selective factors differ")
    entries=map(readlines(f)[2:end]) do line
        row,d,q,name,_,_,_,m,_=split(line,'\t')
        (row=String(row),direction=parse(Int,d),index=parse(Int,q),name=String(name),factor=parse(Int,m)//1000)
    end
    length(entries)==26 && getproperty.(entries,:index)==collect(1:26) || error("Missing odd weights")
    entries
end

struct CleanNumerators{H,L} <: AbstractArray{Interval{BigFloat},3}
    hi::H
    lo::L
end
Base.size(A::CleanNumerators) = size(A.hi)
function Base.getindex(A::CleanNumerators,i::Int,j::Int,k::Int)
    # hi+lo has at most 64 significant bits by construction.
    precision(BigFloat)>=128 || error("Insufficient precision to recover the exact coefficient numerator")
    return interval(BigFloat,BigFloat(A.hi[i,j,k])+BigFloat(A.lo[i,j,k]))
end

function clean_fixed_sources(odd,qo,qe,design)
    design["odd_names"]==odd.names==getproperty.(qo.functionals,:name) || error("Odd sources are out of order")
    design["even_names"]==getproperty.(qe.functionals,:name) || error("Even sources are out of order")
    factors = design["factor_numerators"]
    fixed(A) = begin
        hi,lo = design[A*"_hi"],design[A*"_lo"]
        size(hi)==size(lo)==size(odd.F) || error("Incorrect odd coefficient dimensions")
        raw=getproperty(odd,Symbol(A))
        for j in axes(raw,3), i in CartesianIndices((axes(raw,1),axes(raw,2)))
            x=Float64(inf(raw[i,j])); h=x*factors[j]
            h==hi[i,j] && fma(x,Float64(factors[j]),-h)==lo[i,j] ||
                error("Embedded coefficient differs for $A/$j/$i")
        end
        CleanNumerators(hi,lo)
    end
    odd=merge(odd,(F=fixed("F"),G=fixed("G"),denominators=design["denominators"]))
    for (j,f) in enumerate(qo.functionals)
        scale=factors[j]//design["denominators"][j]
        isequal_interval(interval(Float64,f.nu)*interval(Float64,scale)^2,design["odd_nu"][j]) ||
            error("Wrong fixed odd penalty $j")
    end
    all(isequal_interval.(getproperty.(qe.functionals,:nu),design["even_nu"])) || error("Wrong fixed even penalties")
    weighted(q,k) = merge(q,(functionals=[merge(f,(nu=ν,)) for (f,ν) in zip(q.functionals,design[k])],))
    return odd, weighted(qo,"odd_nu"), weighted(qe,"even_nu")
end

function clean_load_design(folder)
    path=joinpath(folder,"Coercivity",CLEAN_INPUT_NAME)
    d=load(path)
    d["schema"]=="clean-exact-spline-coefficients-v1" && d["factors_sha256"]==CLEAN_WEIGHTS_SHA ||
        error("Incompatible clean input data")
    weights=clean_design(folder)
    d["odd_names"]==getproperty.(weights,:name) &&
        d["denominators"]==denominator.(getproperty.(weights,:factor)) || error("Incorrect coefficient sources or denominators")
    d["base_hashes"]==Dict(n=>clean_hash(joinpath(folder,n)) for n in CLEAN_BASE_FILES) || error("Underlying fixed data changed")
    # The numerators are metadata for validating the encoded fixed coefficients.
    # The numerical kernel uses only their stored limbs and exact denominators.
    d["factor_numerators"] = numerator.(getproperty.(weights,:factor))
    length(d["even_mu"])==length(d["even_nu"])==32 && length(d["odd_nu"])==26 || error("Incomplete fixed weights")
    d
end
