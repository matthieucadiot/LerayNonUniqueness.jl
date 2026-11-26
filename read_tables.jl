# Read the original aggregate certificates of the clean Della verification.
using SHA, TOML, JLD2, IntervalArithmetic, RadiiPolynomial
length(ARGS) == 2 || error("Usage: julia --project=PROOF/residual read_tables.jl EVIDENCE PROOF")
const EVIDENCE, PROOF = abspath.(ARGS)
const RUN = "2d9cec66690641fa88a3aac70215c5cb"
const RECEIPT_SHA = "42fc3d336755a2c10b9d1338ba298e66bc86b7b1cbe8a3be32451b2703292104"
const COERCIVITY = joinpath(EVIDENCE, "results", "coercivity")
const PROP = joinpath(COERCIVITY, "proposition")
verified_file(path, hash) = bytes2hex(open(sha256, path)) == hash ? path : error("SHA-256 mismatch for $(basename(path))")
const RECEIPT = TOML.parsefile(verified_file(joinpath(EVIDENCE,"verification.toml"), RECEIPT_SHA))
@assert RECEIPT["run_id"] == RUN && RECEIPT["all_checks_passed"] === true

# Parse paper decimal targets exactly, independently of binary floating-point.
function decimal(s)
    digits, exp = split(lowercase(s) * (occursin('e',lowercase(s)) ? "" : "e0"), 'e')
    places = occursin('.', digits) ? length(last(split(digits,'.'))) : 0
    parse(BigInt, replace(digits,"."=>"")) * (big(10)//1)^(parse(Int,exp)-places)
end
exact(x) = Rational{BigInt}(x)
valid(x) = isguaranteed(x) && isfinite(inf(x)) && isfinite(sup(x))
encloses(x, d) = valid(x) && exact(inf(x)) <= decimal(d) <= exact(sup(x))

# Outward decimal rounding for the public display (comparisons use full precision).
function disp(x, up=true; digits=9)
    iszero(x) && return "0"
    q = exact(x)
    exponent = floor(Int,log10(abs(Float64(q)))) - digits + 1
    scaled = q / (big(10)//1)^exponent
    number = up ? ceil(BigInt,scaled) : floor(BigInt,scaled)
    s = string(abs(number))
    power = length(s) - 1 + exponent
    sign = number < 0 ? "-" : ""
    if power < -3 || power > 8
        return sign * first(s,1) * "." * s[2:end] * "e" * string(power)
    elseif exponent >= 0
        return sign * s * repeat("0",exponent)
    elseif length(s) + exponent > 0
        return sign * s[1:end+exponent] * "." * s[end+exponent+1:end]
    end
    sign * "0." * repeat("0",-length(s)-exponent) * s
end
row(args...) = println("| ", join(args, " | "), " |")
function table(title, headers...)
    println("\n### ",title,"\n")
    row(headers...)
    row(fill("---",length(headers))...)
end
function upper(name,x,target)
    @assert valid(x) && exact(sup(x)) <= decimal(target) name
    row(name, "≤ " * disp(sup(x)), "≤ " * target)
end
function between(name,x,y,low,high)
    @assert valid(x) && valid(y) && decimal(low) <= exact(inf(x)) <= exact(sup(y)) <= decimal(high) name
    row(name, "["*disp(inf(x),false)*", "*disp(sup(y))*"]", "["*low*", "*high*"]")
end

function profile_tables()
    hash = RECEIPT["upstream"]["residual_sha256"]
    path = verified_file(joinpath(EVIDENCE,"results","residual","Residual_Table1.jld2"),hash)
    jldopen(path) do f
        @assert f["run_id"] == RUN && f["certified"] && f["columns"] == 1:804 && f["rows"] == 1:575
        r, t = f["residual"], f["table1"]
        c = t.profile_constants
        @assert t.certified && c.certified
        for (field,x) in (("odd_residual_upper_exact",r.ℛ),("even_residual_upper_exact",r.ℛ̃))
            @assert exact(sup(x)) == parse(Rational{BigInt},RECEIPT["upstream"][field])
        end
        table("Table 1. Quantitative bounds for Proposition 2.1", "Quantity","Certified bound","Paper")
        upper("Nonlinear residual in Z′",r.ℛ,"6e-12")
        upper("Eigenpair residual in Z′",r.ℛ̃,"1.5e-12")
        upper("Stream residual in Z",c.RΨZ,"1.251e-9")
        between("Even potential gradient in L²",c.gΨ̃_min,c.gΨ̃,"0.88","0.89")
        upper("rΩ* in L∞",r.sΩ,"13")
        upper("Weighted even profile in L²",c.RΩ̃,"22.55")
        upper("Drift field in L∞",t.drift,"22.45")
        upper("Joint L² norm of the two test fields",t.test,"14.4")

        table2()

        table("Table 3. Bounds on residuals and coefficients","Quantity","Certified bound","Paper")
        upper("rΩ* in L∞",r.sΩ,"13")
        upper("rΩ̃* in L∞",r.sΩ̃,"0.025")
        upper("rΩ* in L∞ outside B₉₀₀",r.sΩfar,"5e-4")
        upper("rΩ̃* in L∞ outside B₉₀₀",r.sΩ̃far,"1e-18")
        upper("Iₒ(Ω*, 900)",r.IoΩ,"0.38")
        upper("Iₒ(Ω̃*, 900)",r.IoΩ̃,"3.2e-5")
        between("Ω̃* in L²",r.nΩ̃,r.nΩ̃,"1","1.01")
        upper("√w times the Ω residual in L²",r.wRΩ,"4.3e-12")
        upper("r times the Ω̃ residual in L²",r.rRΩ̃,"1.25e-12")
        upper("√w times the interior Ψ residual in L²",r.wRΨin,"6.8e-14")
        upper("r times the Ψ̃ residual in L²",r.rRΨ̃,"1.6e-14")
        upper("Band Ψ residual in L²",r.RΨband,"5.6e-13")
        upper("√w times the band Ψ residual in L²",r.wRΨband,"1.25e-9")
    end
end

function table2()
    cert=RECEIPT["table2"]
    @assert TOML.parsefile(verified_file(joinpath(COERCIVITY,"table2_checks.toml"),RECEIPT["table2_sha256"])) == cert
    @assert cert["run_id"] == RUN && cert["all_checks_passed"] && cert["physical_passed"]
    f = verified_file(joinpath(COERCIVITY,"Table2_bounds.jld2"),RECEIPT["table2_data_sha256"])
    jldopen(f) do j
        @assert j["run_id"] == RUN && j["counts"] == cert["coverage"]
        bounds = j["table"]
        table("Table 2. Functional and response sums","Parity and sum","Certified upper","Paper upper")
        for (side,targets) in ((:odd,("19.2","19.86","5e-4","15.25")),(:even,("34","61","1e-3","146")))
            for (key,label,target) in zip((:q,:q2,:qe,:Φq),
                    ("functional norms","squared functional norms","functional × residual norms","response × functional norms"),targets)
                x = bounds[side][key]
                rec = cert[String(side)][String(key)]
                @assert rec["passed"] && valid(x) && sup(x) == rec["upper"] && 0 <= exact(sup(x)) <= decimal(target)
                row("$side $label",disp(sup(x)),target)
            end
        end
    end
end

const CELLS = [
    ("odd_T0U0","3.76","12.5455",9), ("odd_T0U1","0.504","1.7237",1),
    ("odd_T1U0","4.96","16.6307",11), ("odd_T1U1","0.489","1.6562",3),
    ("odd_T2U0","0.362","1.3053",1), ("odd_T2U1","0.235","0.854",1),
    ("even_T0U0","3.71","12.3774608652",12), ("even_T0U1","0.494","1.7237277896",0),
    ("even_T1U0","4.72","15.7151282883",9), ("even_T1U1","0.477","1.623936892",2),
    ("even_T2U0","0.285","1.0866964506",0), ("even_T2U1","0.184","0.6501580511",0),
    ("scalar_S0V0","0.706","11.2020998964",6), ("scalar_S1V0","0.706","11.2110651424",3),
    ("scalar_S2V1","0.00238","0.0378391209",0)]
function cell_table()
    p = RECEIPT["proposition"]
    jldopen(verified_file(joinpath(PROP,"Bounds.jld2"),p["bounds_sha256"])) do j
        @assert j["run_id"] == RUN && j["support_bounds"] && j["confinement"] && j["finite_farfield"] && j["farfield"].ok
    end
    jldopen(verified_file(joinpath(PROP,"Boxes.jld2"),p["boxes_sha256"])) do j
        @assert j["run_id"] == RUN && length(j["boxes"].columns) == p["box_columns"] == 2482
    end
    cert = RECEIPT["table4"]
    @assert TOML.parsefile(verified_file(joinpath(PROP,"table4_checks.toml"),RECEIPT["table4_sha256"])) == cert
    @assert cert["run_id"] == RUN && cert["all_checks_passed"] && length(cert["rows"]) == length(CELLS)
    manifest = verified_file(joinpath(PROOF,"CLEAN_SHA256SUMS"),RECEIPT["package_sha256"])
    input_path = "data/Coercivity/Proposition_float64.jld2"
    entry = only(filter(line -> length(line) > 66 && line[67:end] == input_path, readlines(manifest)))
    inputs = JLD2.load(verified_file(joinpath(PROOF,input_path),entry[1:64]),"rows")
    @assert length(inputs) == length(CELLS)
    table("Table 4. Cell construction and matrix certificates","Row","Wᵢ","Λ*","nᵢ","Verified gap ≥","Norm ≤","Least Gershgorin margin ≥")
    for (i,(id,w,target,n)) in enumerate(CELLS)
        input, r = inputs[i], cert["rows"][id]
        @assert input.id == id && input.n == n && encloses(input.W,w) && encloses(input.Λstar,target)
        @assert r["verified"] && r["selected"] == n
        path = verified_file(joinpath(PROP,"matrices",id*".jld2"),r["result_sha256"])
        jldopen(path) do j
            @assert j["row_id"] == id && j["run_id"] == RUN
            c, g, margins = j["certificate"], j["gap"], j["margins"]
            @assert all(valid,(c.normN,c.epsilon,c.tau,g)) && all(valid,margins)
            margin = minimum(inf.(margins))
            @assert length(margins) == r["directions"]
            @assert margin == r["minimum_margin_lower"] && margin >= 0
            @assert sup(c.normN) == r["norm_upper"] && 0 <= sup(c.normN) <= 30
            @assert sup(c.tau) == r["orthogonality_upper"] && 0 <= exact(sup(c.tau)) <= decimal("1e-12")
            @assert sup(c.epsilon) == r["loss_upper"] && 0 <= exact(sup(c.epsilon)) <= decimal("9e-11")
            @assert inf(g) == r["complementary_gap_lower"] == p["spectral_gaps_lower"][i]
            @assert exact(inf(g)) >= decimal(target)
            row(id,w,target,n,disp(inf(g),false;digits=15),disp(sup(c.normN)),disp(margin,false))
        end
    end
end

const WEIGHTS = [
    ("T0",("0.5","1","15","17"),19), ("T1",("15","17","29","31"),21),
    ("T2",("29","31","42","45"),6), ("S0",("0.1","0.2","15","17"),18),
    ("S2",("29","31","152","160"),9), ("U0",("-7","-5","5","7"),18),
    ("U1",("5","7","20","23"),8), ("V0",("-40","-35","35","40"),84),
    ("V1",("-200","-180","180","200"),24)]
function spectral_table()
    path = verified_file(joinpath(EVIDENCE,"results","eigenfunctions","Eigenfunctions.jld2"),RECEIPT["upstream"]["eigenfunctions_sha256"])
    jldopen(path) do j
        @assert j["certified"] && j["run_id"] == RUN
        modes, constants = j["certificates"], j["constants"]
        @assert length(modes) == 207 && length(constants) == 9
        table("Table 5. Spectral weights and verified modes","Weight","(a, x₂, x₃, b)","Highest verified index","Largest eigenvalue error ≤","Largest error in orders 0–2 ≤")
        for (name,knots,m) in WEIGHTS
            subset = filter(c -> c.name == name, modes)
            @assert sort([c.m for c in subset]) == 1:m
            @assert valid(constants[name]) && inf(constants[name]) > 0
            for c in subset
                @assert c.zeros == c.m && c.N == 200 && c.prec == 128
                @assert length(c.nodes) == 6 && length(c.errors) == 3
                @assert all(encloses(c.nodes[i],d) for (i,d) in zip((1,3,4,6),knots))
                @assert all(valid,(c.r,c.Z,c.norm2,c.theta[2],c.errors...))
                @assert inf(c.norm2) > 0 && inf(c.theta[2]) > 0
                @assert exact(sup(c.r)) <= decimal("3.3e-31") && exact(sup(c.Z)) <= decimal("4.1e-12")
                @assert all(e->exact(sup(e)) <= decimal("1.1e-24"),c.errors)
            end
            row(name,"("*join(knots,", ")*")",m,disp(maximum(sup(c.r) for c in subset)),disp(maximum(sup(e) for c in subset for e in c.errors)))
        end
    end
end

# Restore saved 128-bit BigFloat endpoints under their original precision on load.
setprecision(BigFloat,128) do
    profile_tables()
    cell_table()
    spectral_table()
end
