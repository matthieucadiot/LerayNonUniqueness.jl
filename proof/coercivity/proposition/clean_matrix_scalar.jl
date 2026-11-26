# The signed scalar field and its square are independently enclosed by Taylor models.
const field_order = 8
short(x) = functional_order(x,field_order)

function spline_piece(profile,t,lo,hi)
    k=knot_interval(t,lo)
    # Exact BigFloat endpoint enclosures for the local B-splines, then an
    # outward cast. The interval contractions enclose every input coefficient.
    basis=proof_float64(local_basis(t,profile.p,k,lo,hi))
    return (;k,basis,half=proof_float64((hi-lo)/I"2"))
end

function xi_moments(basis,coordinate,physical_range,indices,weight,axis)
    modes=Dict(i=>functional_mode(basis,i,coordinate,physical_range;N=24).exact.f for i in indices)
    bump=functional_bump(coordinate,physical_range,basis.nodes[[1,3,4,6]]).f
    return gram_moments(bump*weight,modes,axis,field_order)
end

function axial_data(inputs,row,a,b,c,d,k_modes)
    f=profile_jets(inputs.profile,(a,b,c,d),12,inputs.far)
    B=f.B
    z2=(B^4-I"1")/I"12"
    c2=I"3"*(B^4-I"1")*B^I"-6"
    c1=I"9"*(B^4+I"1")*B^I"-7"
    vert=spline_piece(inputs.profile,inputs.profile.tv,c,d)
    z=gram_order(f.z,24)
    H=xi_moments(inputs.coercivity.bases[row.z_basis],z,strip_z_range(z,c,d),k_modes,gram_order(f.dz,24),2)
    pairs=column_pairs(H)
    return (;vert,H,pairs,Ht=transpose(column_table(H,pairs,field_order)),hn=[H.norms[p] for p in pairs],
            axis=isthinzero(c),c2=short(c2),c1=short(c1),zfac=short(I"6"*z2*B^I"-3"),
            derivative2=short(-I"6"*B^I"-3"+I"18"*B^I"-7"+c1),
            derivative1=short(-I"27"*B^I"-4"-I"63"*B^I"-8"))
end

function radial_data(inputs,row,a,b,c,d,m_modes)
    f=profile_jets(inputs.profile,(a,b,c,d),12,inputs.far)
    r2=f.r2; A=f.A
    rrange=functional_range(r2)
    inf(rrange)>0 && sup(rrange)<4*inf(rrange) || error("need radial subdivision")
    R=xi_moments(inputs.coercivity.bases[row.r_basis],gram_order(f.r,24),f.Ir,
          m_modes,gram_order(r2,24)*gram_order(A,24)^3/I"6",1)
    return (;rad=spline_piece(inputs.profile,inputs.profile.tu,a,b),R,
            c2=short(I"3"*(A^4-I"1")*A^I"-6"),c1=short(I"3"*(I"5"*A^4+I"3")*A^I"-7"),
            add=short(I"0.035"*r2^I"-1"-row.η))
end

# The six derivative coefficients are those of coefficients_from_jets in
# proposition/Profile.jl. Differentiate the full local spline before reducing
# the Taylor order, and use the production v=s² substitution on the axis.
function scalar_field(profile,r,z)
    p=profile.p
    coeffs=profile.Ψ[r.rad.k-p:r.rad.k,z.vert.k-p:z.vert.k]
    M=r.rad.basis'*coeffs*z.vert.basis
    P=zeros(Interval{Float64},Taylor(p)⊗Taylor(p))
    for i in 0:p, j in 0:p
        P[(i,j)]=M[i+1,j+1]
    end
    function derivative(i,j)
        raw=Model((one(first(coefficients(P)))/(r.rad.half^interval(i)*z.vert.half^interval(j)))*differentiate(P,(i,j)),field_order)
        q=short(raw)
        return z.axis ? short(substitute_axis(q)) : q
    end
    p02=derivative(0,2); p01=derivative(0,1)
    lap=r.c2*derivative(2,0)+r.c1*derivative(1,0)+z.c2*p02+z.c1*p01
    lapv=r.c2*derivative(2,1)+r.c1*derivative(1,1)+z.c2*derivative(0,3)+z.derivative2*p02+z.derivative1*p01
    return lap+z.zfac*lapv+r.add
end

function scalar_column(inputs,row,a,b,V,vertical,m_modes,reference)
    r=radial_data(inputs,row,a,b,V[1],V[2],m_modes)
    pairs=vertical[1].pairs
    S=[fill(I"0",field_order+1,length(pairs)) for _ in 1:2]
    E=[fill(I"0",length(pairs)) for _ in 1:2]
    for (j,z) in enumerate(vertical)
        c,d=V[j:j+1]
        M=scalar_field(inputs.profile,r,z)
        if reference && (j==1 || j==max(1,length(vertical)÷2) || j==length(vertical))
            baseline=only(matrix_field(coefficient_models(inputs.profile,(a,b,c,d),12,inputs.far),row))
            clean_finite(intersect_interval(functional_range(baseline),functional_range(M))) ||
                error("original scalar field disagrees on the audit rectangle")
        end
        jac=z.axis ? proof_float64((b-a)*sqrt(d)/I"4") : proof_float64((b-a)*(d-c)/I"4")
        add_column_models!(S,E,(M,M*M),z.H,pairs,jac,field_order;Ht=z.Ht,h=z.hn)
    end
    rp=column_pairs(r.R)
    return column_matrices(row,rp,pairs,finish_column_models(r.R,rp,S,E,field_order))
end

