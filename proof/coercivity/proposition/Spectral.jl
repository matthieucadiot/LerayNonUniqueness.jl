# Spectral.jl : admissibility, the complementary gaps in Proposition 4.21 (v),
# and a numerical completion of the prescribed negative directions.

mode_eigenvalue(b,m) = m==0 ? I"0" : proof_float64(b.modes[m].lambda)

function admissible_mode(row,bases,a,k)
    z = bases[row.z_basis]
    symmetric = isequal_interval(z.nodes[1],-z.nodes[end]) && isequal_interval(z.nodes[3],-z.nodes[4])
    return !symmetric || (-1)^k==(a==1 ? row.parity : -row.parity)
end

# No functions are changed in a row with n=0. Choose all modes strictly below Lambda_*.
function retained_modes(row,bases,a)
    r = bases[row.r_basis]; z = bases[row.z_basis]
    sigma = Tuple{Int,Int}[]
    for m in 0:length(r.modes), k in 0:length(z.modes)
        admissible_mode(row,bases,a,k) || continue
        λ = mode_eigenvalue(r,m)+mode_eigenvalue(z,k)
        if sup(λ)<inf(row.Λstar)
            push!(sigma,(m,k))
        elseif inf(λ)<sup(row.Λstar)
            error("a mode straddles Lambda_* in $(row.id): ($m,$k)")
        end
    end
    return isempty(sigma) ? zeros(Int,0,2) : [sigma[i][j] for i in eachindex(sigma), j in 1:2]
end

function verify_complement_gap(row,bases)
    r = bases[row.r_basis]; z = bases[row.z_basis]
    # Monotonicity of indexed eigenvalues controls every mode beyond this finite rectangle.
    inf(proof_float64(r.modes[end].lambda))≥sup(row.Λstar) && inf(proof_float64(z.modes[end].lambda))≥sup(row.Λstar) || error("not enough eigenpairs for $(row.id)")
    gap = min(proof_float64(r.modes[end].lambda),proof_float64(z.modes[end].lambda))
    for a in 1:row.s
        sigma = a==1 ? row.sigma1 : row.sigma2
        size(sigma,2)==2 && all(x->0≤x,sigma) || error("invalid retained mode set in $(row.id)")
        selected = Set((sigma[i,1],sigma[i,2]) for i in axes(sigma,1))
        length(selected)==size(sigma,1) || error("a retained mode is repeated in $(row.id)")
        for (m,k) in selected
            m≤length(r.modes) && k≤length(z.modes) && admissible_mode(row,bases,a,k) || error("an inadmissible retained mode occurs in $(row.id)")
        end
        for m in 0:length(r.modes), k in 0:length(z.modes)
            admissible_mode(row,bases,a,k) && !((m,k) in selected) || continue
            λ = mode_eigenvalue(r,m)+mode_eigenvalue(z,k)
            inf(λ)≥sup(row.Λstar) || error("the complementary mode ($m,$k), component $a, violates Lambda_* in $(row.id)")
            gap = min(gap,λ)
        end
    end
    inf(row.d*row.Λstar-row.W)>0 || error("d Lambda_* - W is not positive in $(row.id)")
    return gap
end

function verify_spectral_gaps(data)
    return [verify_complement_gap(row,data.bases) for row in data.rows]
end

# Preserve every supplied xi exactly. Only the remaining columns are chosen numerically.
# The caller certifies I-XX' and X'NX by interval arithmetic before accepting this X.
function complete_coercivity_design(row,N)
    K = size(row.sigma1,1)+size(row.sigma2,1); n = row.n
    size(N)==(K,K) || error("wrong matrix dimension in $(row.id)")
    A = zeros(Interval{Float64},K,n)
    for j in 1:n
        A[:,j] = interval.(Float64,vcat(row.fs[j].xi1,row.fs[j].xi2))
    end
    Q = Matrix(qr(hcat(Float64.(mid.(A)),Matrix{Float64}(I,K,K))).Q)
    Z = Q[:,n+1:K]
    M = Float64.(mid.(N))
    H = Z'*((M+M')/2)*Z
    Z = Z*eigen(Symmetric(H)).vectors
    # Reorthogonalize only the numerical complement; every prescribed xi stays exact.
    Z = Matrix(qr(Z).Q)
    X = hcat(A,interval.(Float64,Z))
    ν = interval.(Float64,[f.nu for f in row.fs])
    return (X=X,L=collect(1:n),nu=ν)
end
