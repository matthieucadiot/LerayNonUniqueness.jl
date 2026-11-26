# Offline creation of the fixed inputs. No original input is overwritten.
include("Response_setup.jl")

function split_numerators(A,weights)
    hi=Array{Float64}(undef,size(A));lo=similar(hi)
    for j in axes(A,3)
        m=numerator(weights[j].factor)
        for i in CartesianIndices((axes(A,1),axes(A,2)))
            coefficient=A[i,j]
            isguaranteed(coefficient) && isthin(coefficient) || error("Nonthin response coefficient")
            x=Float64(inf(coefficient))
            BigFloat(x)==BigFloat(inf(coefficient)) || error("Coefficient is not an exact binary64 number")
            h=x*m
            l=fma(x,Float64(m),-h)
            BigFloat(h)+BigFloat(l)==BigFloat(x)*m || error("Invalid exact coefficient numerator")
            hi[i,j]=h;lo[i,j]=l
        end
    end
    return hi,lo
end

function prepare_fixed_inputs(folder)
    setprecision(BigFloat,128)
    w=clean_design(folder)
    odd=load(joinpath(folder,"responses_odd.jld2"))
    even=load(joinpath(folder,"responses_even.jld2"))
    qo=load(joinpath(folder,"functionals_odd.jld2"));qe=load(joinpath(folder,"functionals_even.jld2"))
    getproperty.(w,:name)==odd["names"]==getproperty.(qo["functionals"],:name) || error("Odd functional order differs")
    even["names"]==getproperty.(qe["functionals"],:name) || error("Even functional order differs")
    Fh,Fl=split_numerators(odd["F"],w);Gh,Gl=split_numerators(odd["G"],w)
    nu=[interval(Float64,f.nu)*interval(Float64,x.factor)^2 for (f,x) in zip(qo["functionals"],w)]
    path=joinpath(folder,"Coercivity",CLEAN_INPUT_NAME)
    temp,io=mktemp(dirname(path));close(io)
    try
        jldsave(temp,true;schema="clean-exact-spline-coefficients-v1",F_hi=Fh,F_lo=Fl,G_hi=Gh,G_lo=Gl,
                denominators=denominator.(getproperty.(w,:factor)),odd_nu=nu,even_nu=qe["nu"],
                even_mu=normalized_even_response_mu(vec(even["mu"])), odd_names=odd["names"],even_names=even["names"],
                factors_sha256=CLEAN_WEIGHTS_SHA,base_hashes=Dict(n=>clean_hash(joinpath(folder,n)) for n in CLEAN_BASE_FILES))
        Base.Filesystem.rename(temp,path)
    finally
        isfile(temp) && rm(temp)
    end
    println("Fixed clean response coefficients written to ",path)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    prepare_fixed_inputs(PROOF_DATA)
end
