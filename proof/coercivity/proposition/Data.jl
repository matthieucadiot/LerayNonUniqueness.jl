# Data.jl : the fifteen rows of Table 4 and the normalized modes of Proposition A.3.
# The selected directions xi and the numbers nu are those of the supplied functionals.

function coercivity_rows()
    names = (("T0","U0"),("T0","U1"),("T1","U0"),("T1","U1"),("T2","U0"),("T2","U1"))
    Wo = ("3.76","0.504","4.96","0.489","0.362","0.235")
    Lo = ("12.5455","1.7237","16.6307","1.6562","1.3053","0.854")
    no = (9,1,11,3,1,1)
    We = ("3.71","0.494","4.72","0.477","0.285","0.184")
    Le = ("12.3774608652","1.7237277896","15.7151282883","1.623936892","1.0866964506","0.6501580511")
    ne = (12,0,9,2,0,0)
    rows = NamedTuple[]
    for (block,parity,W,Λ,n) in ((:odd,-1,Wo,Lo,no),(:even,1,We,Le,ne)), i in 1:6
        r,z = names[i]; cell = r*z
        push!(rows,(id=string(block)*"_"*cell,block=block,cell=cell,r_basis=r,z_basis=z,parity=parity,
                    s=2,W=parse(Interval{Float64},W[i]),Λstar=parse(Interval{Float64},Λ[i]),n=n[i],
                    d=I"0.899999",η=I"0.01"))
    end
    for (cell,r,z,W,Λ,n) in (("S0V0","S0","V0","0.706","11.2020998964",6),
                              ("S1V0","T1","V0","0.706","11.2110651424",3),
                              ("S2V1","S2","V1","0.00238","0.0378391209",0))
        push!(rows,(id="scalar_"*cell,block=:scalar,cell=cell,r_basis=r,z_basis=z,parity=1,
                    s=1,W=parse(Interval{Float64},W),Λstar=parse(Interval{Float64},Λ),n=n,
                    d=I"0.189",η=I"1e-8"))
    end
    return rows
end

function spectral_weights()
    return (("T0",("0.5","1","15","17"),true,19), ("T1",("15","17","29","31"),true,21),
            ("T2",("29","31","42","45"),true,6), ("S0",("0.1","0.2","15","17"),true,18),
            ("S2",("29","31","152","160"),true,9), ("U0",("-7","-5","5","7"),false,18),
            ("U1",("5","7","20","23"),false,8), ("V0",("-40","-35","35","40"),false,84),
            ("V1",("-200","-180","180","200"),false,24))
end

# Exact knots and the chosen mode parameters stay in the serialized precision.
# Matrix entries and bounds use outward Float64 enclosures throughout.
function load_proof_inputs(file=joinpath(PROOF_DATA,"Proposition_4_21_inputs.jld2"))
    setrounding(BigFloat,RoundNearest)
    r = load_coercivity_proposition(file)
    r["prec"]==128 || error("the source proof inputs require 128-bit interval data")
    rows = [merge(row,(W=proof_float64(row.W),Λstar=proof_float64(row.Λstar),
                     d=proof_float64(row.d),η=proof_float64(row.η))) for row in r["rows"]]
    bases = r["bases"]
    odd = (functionals=vcat([row.fs for row in rows if row.block==:odd]...),eigenfunctions=bases,odd=true,published_mode_bounds=true)
    even = (functionals=vcat([row.fs for row in rows if row.block!=:odd]...),eigenfunctions=bases,odd=false,published_mode_bounds=true)
    return (profile=r["profile"],coercivity=(odd=odd,even=even,bases=bases,rows=rows),far=proof_float64(r["far"]))
end
