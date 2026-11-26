# Response_moments.jl : the same scalar moment integrals as Taylor.jl, contracted before forming a polynomial square.
# If A_ij and B_kl are polynomial coefficients, ∫P_A P_B ρ μ=Σ A_ij (M_ρ B M_μᵀ)_ij,
# (M_ρ)_ik=m_(i+k), (M_μ)_jl=n_(j+l). Coefficient signs and all interval roundoff are retained.

function response_moment_matrix(ρ,n,m)
    length(ρ.m)>n+m || error("not enough moments for the polynomial contraction")
    return [ρ.m[i+k+1] for i in 0:n,k in 0:m]
end

function response_coefficient_matrix(Π)
    n,m=order(space(Π))
    return [Π[(i,j)] for i in 0:n,j in 0:m]
end

function response_contraction(A,B,ρ,μ)
    n1,m1=size(A);n2,m2=size(B)
    Mρ=response_moment_matrix(ρ,n1-1,n2-1)
    Mμ=response_moment_matrix(μ,m1-1,m2-1)
    return sum(A.*(Mρ*B*transpose(Mμ)))
end

function response_unit_moments(d,c)
    return (m=[iseven(k) ? I"2"*one(c)/interval(BigFloat,k+1) : zero(c) for k in 0:d],)
end

function response_integral(F::Model,ρ,μ)
    A=response_coefficient_matrix(F.Π)
    n,m=size(A)
    value=response_contraction(A,A,ρ,μ)
    unitρ=response_unit_moments(2*(n-1),F.Π[(0,0)])
    unitμ=response_unit_moments(2*(m-1),F.Π[(0,0)])
    unit=intersect_interval(response_contraction(A,A,unitρ,unitμ),I"[0,inf]")
    ε=I"2"*F.ϱ*norm(F.Π,1)+F.ϱ*F.ϱ
    δ=ρ.ϱ*μ.Π1+μ.ϱ*ρ.Π1+ρ.ϱ*μ.ϱ
    err=δ*unit+ε*(ρ.m[1]*μ.m[1]+I"4"*δ)
    return intersect_interval(value+I"[-1,1]"*err,I"[0,inf]")
end

function response_pairing(F::Model,G::Model,ρ,μ)
    A=response_coefficient_matrix(F.Π);B=response_coefficient_matrix(G.Π)
    value=response_contraction(A,B,ρ,μ)
    ε=F.ϱ*norm(G.Π,1)+G.ϱ*norm(F.Π,1)+F.ϱ*G.ϱ
    δ=ρ.ϱ*μ.Π1+μ.ϱ*ρ.Π1+ρ.ϱ*μ.ϱ
    # ||Π_F Π_G||₁≤||Π_F||₁||Π_G||₁ avoids constructing the convolution only for its error bound.
    err=I"4"*δ*norm(F.Π,1)*norm(G.Π,1)+ε*(ρ.m[1]*μ.m[1]+I"4"*δ)
    return value+I"[-1,1]"*err
end
