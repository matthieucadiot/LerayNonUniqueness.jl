# Accept only numerical results produced in the given clean run.
using SHA, TOML, Dates, IntervalArithmetic, JLD2
include("clean_closure.jl")
module Residual
include(joinpath(@__DIR__,"..","residual","Setup.jl"))
end
module Proposition
include(joinpath(@__DIR__,"..","coercivity","proposition","Setup.jl"))
include(joinpath(@__DIR__,"..","coercivity","proposition","Files.jl"))
end
module Prerequisites
include(joinpath(@__DIR__,"..","coercivity","Certificates.jl"))
end
digest(path) = bytes2hex(open(sha256,path))

function check_prerequisites()
    profilefile = ENV["RESIDUAL_CERTIFICATE"]
    Prerequisites.check_shared_sources()
    eig = Prerequisites.load_spectral_certificate()
    cert = JLD2.load(profilefile)
    cert["run_id"]==ENV["CLEAN_RUN_ID"] && cert["signature"]==Residual.residual_signature() ||
        error("Residual certificate run or input mismatch")
    cert["columns"]==1:804 && cert["rows"]==1:575 || error("Incomplete residual grid")
    Residual.verify_bounds(cert["residual"])
    Residual.verify_table1(cert["table1"])
    residual = cert["residual"]
    sup(residual.ℛ) <= big(6)//big(10)^12 && sup(residual.ℛ̃) <= big(15)//big(10)^13 ||
        error("The fixed-point residual inputs exceed the proof targets")
    Dict("eigenfunctions_sha256"=>eig.sha256,"residual_sha256"=>digest(profilefile),
         "nonconstant_modes"=>207,"constant_modes"=>9,"profile_columns"=>804,"profile_rows"=>575,
         "odd_residual_upper_exact"=>string(Rational{BigInt}(sup(residual.ℛ))),
         "even_residual_upper_exact"=>string(Rational{BigInt}(sup(residual.ℛ̃))))
end

function check_coefficient_evidence()
    P = Proposition
    inputs=P.load_proof_inputs()
    U=P.box_partition(inputs.profile.U,2482)
    V=P.box_partition(inputs.profile.V)
    file=joinpath(P.OUTPUT,"Boxes.jld2")
    boxes=P.read_proof(file)["boxes"]
    columns=joinpath(P.OUTPUT,"box_columns")
    # assemble_boxes checks each of the 2482 stored partitions, heights and inequalities.
    checked=P.assemble_boxes(columns,U,V)
    length(boxes.columns)==length(checked.columns) && all(i->begin
        a,b=checked.columns[i],boxes.columns[i]
        a.j1==b.j1 && a.j2==b.j2 &&
            all(k->isequal_interval(getproperty(a,k),getproperty(b,k)),(:L1,:L2,:h1,:h2))
    end, eachindex(checked.columns)) || error("Aggregate boxes differ from recomputed heights")
    length(boxes.r)==length(checked.r) && all(isequal_interval.(checked.r,boxes.r)) || error("Box radii differ")
    length(boxes.columns)==2482 && length(boxes.U)==length(U) && length(boxes.V)==length(V) &&
        all(isequal_interval.(boxes.U,U)) && all(isequal_interval.(boxes.V,V)) || error("Box grid mismatch")
    C=P.read_proof(joinpath(P.OUTPUT,"Bounds.jld2"))
    all(k->C[k]===true,("finite_farfield","confinement","support_bounds")) && C["farfield"].ok === true ||
        error("Missing coefficient acceptance")
    p=inputs.profile
    exterior=P.exterior_farfield(sqrt(P.U0(interval(BigFloat,p.L))))
    saved=C["farfield"]
    exterior.ok === true && isequal_interval(exterior.radius,saved.radius) &&
        keys(exterior.bounds)==keys(saved.bounds) &&
        all(k->isequal_interval(exterior.bounds[k],saved.bounds[k]),keys(exterior.bounds)) ||
        error("Analytic exterior failed or differs from its saved certificate")
    x,y=P.bounds_grid(inputs,boxes)
    for (kind,nu,nv) in (("farfield",length(p.U)-1,length(p.V)-1),("bounds",length(x)-1,length(y)-1))
        for i=1:nu
            col=P.read_proof(joinpath(P.OUTPUT,kind*"_columns","column_$i.jld2"))
            col["verified"]===true && col["column"]==i && col["rows"]==1:nv ||
                error("Incomplete $kind coverage at $i")
        end
    end
    gaps=P.verify_spectral_gaps(inputs.coercivity)
    length(gaps)==15 || error("Missing complementary gaps")
    Dict("box_columns"=>length(boxes.columns),"support_columns"=>length(x)-1,"support_rows"=>length(y)-1,
         "farfield_columns"=>length(p.U)-1,"farfield_rows"=>length(p.V)-1,
         "spectral_gaps_lower"=>Float64.(inf.(gaps)),
         "boxes_sha256"=>digest(file),"bounds_sha256"=>digest(joinpath(P.OUTPUT,"Bounds.jld2")))
end

function final()
    run=realpath(ENV["CLEAN_RUN_DIR"])
    result=realpath(ENV["COERCIVITY_RESULTS"])
    result==joinpath(run,"results","coercivity") || error("Results must belong to this run")
    upstream=check_prerequisites()
    coeff=check_coefficient_evidence()
    tf=joinpath(result,"table2_checks.toml"); mf=joinpath(result,"proposition","table4_checks.toml")
    t=TOML.parsefile(tf); m=TOML.parsefile(mf)
    t["run_id"]==m["run_id"]==ENV["CLEAN_RUN_ID"] || error("Component reports belong to different runs")
    t["all_checks_passed"]===true && t["physical_passed"]===true && t["coverage"]["columns"]==67682 &&
        t["coverage"]["parent_cells"]==57725138 && t["coverage"]["selected_columns"]==4289 || error("Incomplete response coverage")
    all(k->t[k]==upstream[k],("eigenfunctions_sha256","residual_sha256")) || error("Response prerequisites differ")
    t["table2_data_sha256"]==digest(joinpath(result,"Table2_bounds.jld2")) &&
        t["fixed_coefficients_sha256"]==digest(joinpath(ENV["PROOF_DATA"],"Coercivity","Fixed_weighted_responses.jld2")) ||
        error("Response results or fixed coefficient file differs from the assembled certificate")
    targets=(odd=(q=96//5,q2=993//50,qe=1//2000,Φq=61//4),even=(q=34//1,q2=61//1,qe=1//1000,Φq=146//1))
    for side in keys(targets)
        for (key,target) in pairs(targets[side])
            rec=t[string(side)][string(key)]
            rec["passed"]===true && 0 <= rec["upper"] <= target || error("Failed Table 2 $side $key")
        end
        norms=t["physical"][string(side)]
        length(norms["individual_B_squared_upper"])==(side==:odd ? 26 : 32) &&
            0<=norms["sum_B_upper"]<=targets[side].q &&
            0<=norms["sum_B_squared_upper"]<=targets[side].q2 || error("Functional sums failed")
    end
    rows=m["rows"]
    expected=Set(r.id for r in Proposition.load_proof_inputs().coercivity.rows)
    m["all_checks_passed"]===true && m["matrix_count"]==15 && Set(keys(rows))==expected || error("Incomplete Table 4")
    for (id,row) in rows
        matrixfile=joinpath(result,"proposition","matrices",id*".jld2")
        row["result_sha256"]==digest(matrixfile) && row["verified"]===true || error("Matrix record differs for $id")
        0<=row["norm_upper"]<=30 && 0<=row["loss_upper"]<=big(9)//10^11 &&
            row["minimum_margin_lower"]>=0 || error("Matrix inequalities failed for $id")
    end
    closing=closure()
    length(closing)==7 || error("Missing closing checks")
    report=Dict("schema"=>"leray-clean-certification-v1","all_checks_passed"=>true,"run_id"=>ENV["CLEAN_RUN_ID"],
                "completed_utc"=>string(now(UTC)),"julia"=>string(VERSION),
                "package_sha256"=>digest(joinpath(@__DIR__,"..","CLEAN_SHA256SUMS")),
                "upstream"=>upstream,"proposition"=>coeff,"table2"=>t,"table4"=>m,
                "table2_data_sha256"=>t["table2_data_sha256"],
                "fixed_weights_sha256"=>t["fixed_coefficients_sha256"],
                "table2_sha256"=>digest(tf),"table4_sha256"=>digest(mf),"closure"=>closing)
    file=joinpath(run,"verification.toml")
    tmp,io=mktemp(run)
    try
        TOML.print(io,report;sorted=true); close(io)
        Base.Filesystem.rename(tmp,file)
    finally
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp)
    end
    println("CLEAN CERTIFICATION PASSED ",digest(file),"  ",file)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && final()
