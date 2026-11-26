# Response_data.jl : the response representation (6.24), on the required spline refinements.
# Input Float64 interval coefficients are thin exact binary numbers defining the chosen approximate responses.

# This prescribed binary multiplier selects new auxiliary scalars for the fixed,
# already normalized eigenfunction. It does not change any profile or input file.
const RESPONSE_EVEN_MU_MULTIPLIER = 0.9952979094776627

function normalized_even_response_mu(μ)
    # The product of two Float64 numbers is exact at 256-bit precision. Round it
    # once to the chosen Float64 coefficient, then enclose that coefficient as
    # a thin interval; no assertion about an exact normalization ratio is used.
    return setprecision(BigFloat,256) do
        setrounding(BigFloat,RoundNearest) do
            factor=BigFloat(RESPONSE_EVEN_MU_MULTIPLIER)
            map(μ) do x
                y=Float64(BigFloat(inf(x))*factor,RoundNearest)
                isfinite(y) || error("a normalized even response scalar is not finite")
                interval(Float64,y)
            end
        end
    end
end

function check_jld2_complete(file)
    actual=filesize(file)
    signature=UInt8[0x89,0x48,0x44,0x46,0x0d,0x0a,0x1a,0x0a]
    open(file) do io
        for offset in (0,512,1024,2048,4096)
            actual≥offset+48 || continue
            seek(io,offset)
            read(io,8)==signature || continue
            version=read(io,UInt8); width=read(io,UInt8)
            if version in (2,3) && width==8
                seek(io,offset+28)
                expected=ltoh(read(io,UInt64))
                expected≤actual || error("$file is truncated: $actual bytes present, but its HDF5 header requires $expected bytes")
            end
            return true
        end
        error("$file has no supported HDF5/JLD2 header")
    end
end

function response_integer(x,name)
    y=x isa Interval ? (isthin(x) ? inf(x) : error("$name is not exact")) : x
    y isa Real && isinteger(y) || error("$name is not an integer")
    return Int(y)
end

function validate_response_data(r,odd::Bool; embedded_mu=nothing)
    needed=odd ? ("F","G","u_knots","v_knots","degree","L","d","alpha_F","alpha_G","names") :
                 ("F","G","u_knots","v_knots","degree","L","d","alpha_F","alpha_G","mu","a","p","names")
    all(k->haskey(r,k),needed) || error("the response file needs keys $needed; its keys are $(collect(keys(r)))")
    p=response_integer(r["degree"],"degree");L=response_integer(r["L"],"L");d=response_integer(r["d"],"d")
    αF=response_integer(r["alpha_F"],"alpha_F");αG=response_integer(r["alpha_G"],"alpha_G")
    (p,L,d)==(9,100,16) && (αF,αG)==(odd ? (4,2) : (3,1)) || error("the response parameters differ from (6.24)")
    tu=interval.(BigFloat,vec(r["u_knots"]));tv=interval.(BigFloat,vec(r["v_knots"]))
    n=odd ? 26 : 32
    F=r["F"];G=r["G"]
    size(F)==size(G)==(length(tu)-p-1,length(tv)-p-1,n) || error("the response coefficient dimensions are incompatible")
    for A in (tu,tv,F,G)
        all(x->x isa Interval && isguaranteed(x)&&isthin(x)&&isfinite(inf(x)),A) || error("response knots and coefficients must be finite guaranteed thin intervals")
    end
    for t in (tu,tv)
        issorted(inf.(t)) && all(x->inf(x)==0,t[1:p+1]) && all(x->inf(x)==L,t[end-p:end]) || error("response knots must be sorted and clamped")
        for x in unique(inf.(t))
            (x==0 || x==L || count(y->inf(y)==x,t)≤p-1) || error("response splines must be continuously differentiable")
        end
    end
    μ=odd ? fill(I"0",n) : interval.(Float64,vec(r["mu"]))
    a=odd ? fill(I"0",n,4) : interval.(Float64,r["a"])
    length(μ)==n && size(a)==(n,4) || error("the even scalar/multipole dimensions are incompatible")
    all(x->isguaranteed(x)&&isthin(x)&&isfinite(inf(x)),μ) && all(x->isguaranteed(x)&&isthin(x)&&isfinite(inf(x)),a) || error("the response scalar and multipole coefficients are not exact finite intervals")
    odd || r["p"]==[[1],[4,-5],[8,-28,21],[64,-432,792,-429]] || error("the multipole polynomials differ from (6.23)")
    names=String.(r["names"])
    length(names)==n && length(unique(names))==n || error("the response names are inconsistent")
    # Always start from the original stored scalars. The input dictionary and
    # its arrays are not mutated, so repeated loads do not accumulate scaling.
    if !odd
        reference=normalized_even_response_mu(μ)
        if isnothing(embedded_mu)
            μ=reference
        else
            all(isequal_interval.(reference,embedded_mu)) || error("Embedded even scalars differ from the fixed basis normalization")
            μ=embedded_mu
        end
    end
    length(μ)==n && all(x->isguaranteed(x) && isthin(x) && isfinite(inf(x)),μ) || error("invalid fixed response scalars")
    return (p=p,L=L,d=d,tu=tu,tv=tv,F=F,G=G,μ=μ,a=a,n=n,odd=odd,alpha_F=αF,alpha_G=αG,names=names)
end

function load_response_data(file,odd::Bool)
    check_jld2_complete(file)
    return validate_response_data(load(file),odd)
end

# Thin dyadic grid endpoints: include every spline knot and the cutoff boundary exactly.
function response_partition(profile,odd,even,functionals=())
    U=collect(inf.(profile.U));V=collect(inf.(profile.V))
    append!(U,inf.(odd.tu));append!(U,inf.(even.tu))
    append!(V,inf.(odd.tv));append!(V,inf.(even.tv))
    push!(U,BigFloat(profile.L-profile.d));push!(V,BigFloat(profile.L-profile.d))
    for dataset in functionals, b in values(dataset.eigenfunctions), x in b.nodes
        s=abs(interval(BigFloat,x))
        u=(I"1"+I"12"*s*s)^I"0.25"-I"1"
        # Both outward endpoints are exact dyadic integration boundaries. The tiny cell between them
        # covers an irrational mode junction; functional_models encloses both pieces there.
        for y in (inf(u),sup(u))
            0≤y≤profile.L && push!(b.radial ? U : V,y)
        end
    end
    return (U=interval.(BigFloat,sort(unique(U))),V=interval.(BigFloat,sort(unique(V))))
end

include(joinpath(@__DIR__,"Response_grids.jl"))

# Proposition A.3 supplies the published C² error, which is deliberately used instead of trusting
# smaller per-mode numbers that would require the source spectral certificates to be checked again.
function load_functional_data(file,odd::Bool)
    precision(BigFloat)==128 || error("functional interval data must be read at 128 bits")
    setrounding(BigFloat,RoundNearest)
    check_jld2_complete(file)
    r=load(file)
    all(k->haskey(r,k),("functionals","eigenfunctions","nu")) || error("incomplete functional schema")
    fs=r["functionals"];n=odd ? 26 : 32
    length(fs)==length(r["nu"])==n || error("incorrect number of functionals")
    highest=Dict("T0"=>19,"T1"=>21,"T2"=>6,"S0"=>18,"S2"=>9,"U0"=>18,"U1"=>8,"V0"=>84,"V1"=>24)
    plateaux=Dict("T0"=>(1//2,1,15,17),"T1"=>(15,17,29,31),"T2"=>(29,31,42,45),
                 "S0"=>(1//10,1//5,15,17),"S2"=>(29,31,152,160),"U0"=>(-7,-5,5,7),
                 "U1"=>(5,7,20,23),"V0"=>(-40,-35,35,40),"V1"=>(-200,-180,180,200))
    ε=sup(big(I"11")/(big(I"10")^I"25"));radius=big(I"33")/(big(I"10")^I"32")
    bases=copy(r["eigenfunctions"])
    for (name,b) in bases
        haskey(highest,name) && b.N==200 && length(b.nodes)==6 || error("unknown basis in Proposition A.3: $name")
        b.radial==(name[1] in ('T','S')) && all(in_interval.(collect(plateaux[name]),b.nodes[[1,3,4,6]])) || error("basis $name differs from Table 5")
        all(x->isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x)),b.nodes) && all(diff(inf.(b.nodes)).>0) || error("invalid mode patch endpoints")
        modes=[begin
            m.m==i && i≤highest[name] && isguaranteed(m.lambda_bar) && isthin(m.lambda_bar) && inf(m.lambda_bar)>0 || error("mode $name/$i is outside Proposition A.3")
            length(m.errors)==3 && all(x->isfinite(x)&&x≥0,m.errors) || error("invalid mode error bounds")
            inf(m.phi_a[1])>ε || error("the supplied mode does not have the canonical positive endpoint phase")
            for c in (m.phi_a,m.psi_a,m.psi_b,m.phi_b)
                length(c)==201 && all(x->isguaranteed(x)&&isthin(x)&&isfinite(inf(x)),c) || error("invalid spectral polynomial coefficients")
            end
            λ=interval(BigFloat,m.lambda_bar)+I"[-1,1]"*radius
            merge(m,(lambda=λ,errors=[max(x,ε) for x in m.errors]))
        end for (i,m) in enumerate(b.modes)]
        if b.radial
            X=functional_constant(big(I"1"),Model(Sequence(Taylor(0) ⊗ Taylor(0),[big(I"1")]),24))
            mass=functional_bump_integral(X,b.nodes[[1,3,4,6]],4).Π[(0,0)]
        else
            mass=(b.nodes[end]-b.nodes[1]+b.nodes[4]-b.nodes[3])/I"2"
        end
        c0=mass^I"-0.5"
        isempty_interval(intersect_interval(c0,b.constant)) && error("the constant mode normalization is inconsistent")
        bases[name]=merge(b,(modes=modes,constant=c0))
    end
    for (j,f) in enumerate(fs)
        f.parity==(odd ? -1 : 1) && f.kind in ("vector","scalar") || error("invalid functional parity/type")
        isequal_interval(f.nu,r["nu"][j]) && inf(f.nu)≥0 || error("functional nu ordering or sign differs")
        all(x->isguaranteed(x)&&isthin(x)&&isfinite(inf(x)),vcat(f.xi1,f.xi2,[f.nu])) || error("the functional coefficients must be exact binary numbers")
        size(f.sigma1)==(length(f.xi1),2) && size(f.sigma2)==(length(f.xi2),2) || error("functional mode dimensions differ")
        for (name,spec) in ((f.r_basis,f.r_spec),(f.z_basis,f.z_spec))
            haskey(bases,name) || error("unknown functional basis")
            all(isequal_interval.(spec,bases[name].nodes[[1,3,4,6]])) || error("the functional and mode plateaux differ")
        end
        f.symmetric_cell && begin
            all(k->(-1)^k==f.parity,f.sigma1[:,2]) || error("Xi1 has the wrong parity")
            all(k->(-1)^k==-f.parity,f.sigma2[:,2]) || error("Xi2 has the wrong parity")
        end
    end
    return (functionals=fs,eigenfunctions=bases,odd=odd,published_mode_bounds=true)
end

# The positive-half-domain formula uses either an even plateau or the reflected
# copy of a positive-z cell. The compact q support must be inside the integration square.
function validate_functional_support(f,L)
    inf(f.r_spec[1])>0 && sup(f.r_spec[4])≤inf(L) || error("the functional radial support is outside the integration domain")
    max(sup(abs(f.z_spec[1])),sup(abs(f.z_spec[4])))≤inf(L) || error("the functional vertical support is outside the integration domain")
    if f.symmetric_cell
        all(isequal_interval.(f.z_spec,-reverse(f.z_spec))) || error("the symmetric functional plateau is not even")
    else
        inf(f.z_spec[1])≥0 || error("a reflected functional must be specified on the positive vertical half-domain")
    end
    return true
end

function load_table2_data(folder; clean=clean_load_design(folder))
    profile=load_coercivity_profile(joinpath(folder,"Profile"))
    odd=load_response_data(joinpath(folder,"responses_odd.jld2"),true)
    even_file=joinpath(folder,"responses_even.jld2")
    even=isnothing(clean) ? load_response_data(even_file,false) : validate_response_data(load(even_file),false; embedded_mu=clean["even_mu"])
    qo=load_functional_data(joinpath(folder,"functionals_odd.jld2"),true)
    qe=load_functional_data(joinpath(folder,"functionals_even.jld2"),false)
    if !isnothing(clean)
        odd,qo,qe = clean_fixed_sources(odd,qo,qe,clean)
    end
    L=sqrt(U0(interval(BigFloat,profile.L)))
    for q in (qo,qe),f in q.functionals
        validate_functional_support(f,L)
    end
    odd.names==[f.name for f in qo.functionals] && even.names==[f.name for f in qe.functionals] || error("response and functional orderings differ")
    grid=response_partition(profile,odd,even,(qo,qe))
    grids=response_grids(profile,odd,even,(qo,qe))
    on_grid(g)=merge(profile,(U=g.U,V=g.V,λ̃=interval(Float64,profile.λ̃)))
    return (profile=on_grid(grid),odd_profile=on_grid(grids.odd),even_profile=on_grid(grids.even),
            original_profile=profile,odd=odd,even=even,qo=qo,qe=qe)
end

function response_basis_column(i,data,U)
    a=U[i];b=U[i+1];k=knot_interval(data.tu,a)
    return (B=local_basis(data.tu,data.p,k,a,b),k=k,h=interval(Float64,(b-a)/I"2"))
end
function response_basis_row(j,data,V)
    a=V[j];b=V[j+1];k=knot_interval(data.tv,a)
    return (B=local_basis(data.tv,data.p,k,a,b),k=k,h=interval(Float64,(b-a)/I"2"))
end
function response_spline_jets(j,data,C,R,N)
    denominator = hasproperty(data,:denominators) ? data.denominators[j] : 1
    F=field(piece(view(data.F,:,:,j),C.B,C.k,R.B,R.k,data.p,denominator),C.h,R.h,N)
    G=field(piece(view(data.G,:,:,j),C.B,C.k,R.B,R.k,data.p,denominator),C.h,R.h,N)
    return G,F
end
