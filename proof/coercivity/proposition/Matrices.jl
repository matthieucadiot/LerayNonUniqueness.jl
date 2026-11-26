using IntervalArithmetic, RadiiPolynomial, LinearAlgebra


#################### Finite matrices: (4.48) and (4.50) ####################

function check_matrix(A)
    size(A,1)==size(A,2) && size(A,1)>0 || error("A finite matrix must be square and nonempty")
    all(x -> isguaranteed(x) && isfinite(inf(x)) && isfinite(sup(x)), A) ||
        error("A matrix entry is not a finite guaranteed interval")
    all(isequal_interval(A[p,q],A[q,p]) for p in axes(A,1), q in axes(A,2)) ||
        error("The interval matrix must be symmetric")
    return nothing
end

# Lambda contains the retained eigenvalue sums, as a diagonal matrix or a vector.
function assemble_N(d,Λ,G,G2,Λstar,W)
    d,Λstar,W = interval(Float64,d),interval(Float64,Λstar),interval(Float64,W)
    all(isguaranteed,(d,Λstar,W)) || error("A matrix parameter lost its interval guarantee")
    inf(d)>0 && inf(Λstar)>0 && inf(W)>=0 || error("Invalid diffusion, gap, or potential bound")
    G,G2 = interval.(Float64,G),interval.(Float64,G2)
    Λ = Λ isa AbstractVector ? Matrix(Diagonal(interval.(Float64,Λ))) : interval.(Float64,Λ)
    check_matrix(Λ)
    check_matrix(G)
    check_matrix(G2)
    size(Λ)==size(G)==size(G2) || error("The finite matrices have different sizes")
    all(isthinzero(Λ[p,q]) for p in axes(Λ,1), q in axes(Λ,2) if p!=q) ||
        error("Lambda must be diagonal")
    denominator = d*Λstar-W
    inf(denominator)>0 || error("d*Lambda_star-W must be positive")
    N = d.*Λ+G-(G2-G*G)./denominator
    check_matrix(N)
    return N
end

# For a symmetric real matrix, its maximum absolute row sum bounds its spectral norm.
function matrix_row_bound(A)
    return maximum(sum(abs.(A);dims=2))
end

# X is the stored real candidate defining the functionals; it is never recomputed.
# nu[p] belongs to the column L[p], with no ordering requirement on L.
function certify_matrix(N,X,L,ν,η)
    N,X,ν,η = interval.(Float64,N),interval.(Float64,X),interval.(Float64,ν),interval(Float64,η)
    check_matrix(N)
    K = size(N,1)
    size(X)==(K,K) || error("X and N have different sizes")
    all(x -> isguaranteed(x) && isthin(x) && isfinite(inf(x)), X) ||
        error("X must enclose an exact finite real matrix in thin intervals")
    X = convert.(eltype(N),X)
    all(isthin,X) || error("The working precision cannot represent the stored X exactly")
    length(L)==length(ν) && length(unique(L))==length(L) || error("Invalid functional indices")
    all(l -> l isa Integer && 1<=l<=K, L) || error("A functional index is outside 1:K")
    all(x -> isguaranteed(x) && inf(x)>0 && isfinite(sup(x)), ν) || error("nu must be positive and finite")
    isguaranteed(η) && inf(η)>=0 && isfinite(sup(η)) || error("eta must be nonnegative and finite")

    Id = fill(zero(first(X)),K,K)
    for p in 1:K
        Id[p,p] = one(first(X))
    end
    defect = Id-X*transpose(X)
    # The Frobenius norm encloses an upper bound for the spectral defect tau.
    tau = sqrt(sum(abs.(defect).^2))
    tolerance = interval(Float64,1)/interval(Float64,10)^12
    sup(tau)<=inf(tolerance) || error("The orthogonality defect exceeds 1e-12")

    D = transpose(X)*N*X
    # ||N|| <= ||X^T N X||/(1-tau), since the smallest singular value
    # of X is at least sqrt(1-tau). D is nearly diagonal for the stored X.
    normN = min(matrix_row_bound(N),matrix_row_bound(D)/(one(tau)-tau))
    sup(normN)<=30 || error("The spectral norm bound exceeds 30")
    g = [D[l,l]-sum((abs(D[l,k]) for k in 1:K if k!=l);init=zero(D[l,l])) for l in 1:K]
    for l in 1:K
        p = findfirst(==(l),L)
        if p===nothing
            inf(g[l])>=0 || error("A Gershgorin number outside L is not nonnegative")
        else
            inf(g[l]+ν[p])>=0 || error("A Gershgorin number is smaller than -nu")
        end
    end
    epsilon = interval(3)*tau*normN
    all(isguaranteed,(tau,normN,epsilon)) && all(isguaranteed,g) ||
        error("A matrix certificate lost its interval guarantee")
    sup(epsilon)<=inf(η) || error("The loss 3*tau*normN exceeds eta")
    return (tau=tau,normN=normN,D=D,g=g,epsilon=epsilon)
end
