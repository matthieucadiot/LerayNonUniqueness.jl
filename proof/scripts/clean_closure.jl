# Rational closing inequalities of the selective-weight proof.
require(ok,description) = ok === true || error(description)
function closure()
    Q(a,b=1)=big(a)//big(b)
    γ=Q(1,10); δ=Q(1,10^6); κ=Q(1,2000); ρ=Q(1,10^9)
    require(κ<=γ/100,"residual fraction")
    require(δ<Q(1,100)/(6+2*26+2*Q(1986,100)) && δ<Q(1,100)/(6+2*26+2*61)/(1+(Q(135,100)+27)/Q(88,100)),"inverse delta restrictions")
    z=(1+Q(1525,100))/(γ-κ); z0=(1+57δ)*z+2δ
    h=Q(3,10^7); μ=Q(6,10^6)
    ze=(1+146)/(Q(98,100)*γ)
    require(ze==1500 && Q(175,100)+Q(186,10)*ze<28000,"even inverse constants")
    f=Q(191,10^12)
    require(Q(15,10^13)+Q(1888,10000)*ρ<f,"effective even forcing")
    forcing=f+(Q(16,100)*ρ+μ)*h
    λ=Q(235059685902438,10^15)
    require(Q(23,100)<λ-μ<λ+μ<Q(24,100) && h<1,"eigenvalue and nontriviality")
    results=(odd_inverse_Z=(z,Q(164)),odd_inverse_Z0=(z0,Q(164)),
      odd_image=(164*(Q(6,10^12)+Q(8,100)*ρ^2),ρ),odd_Lipschitz=(2*164*Q(8,100)*ρ,Q(1)),
      even_h_image=(ze*forcing,h),even_mu_image=(28000*forcing,μ),
      even_product_Lipschitz=(h*(Q(16,100)*ρ+2μ)*max(ze/h,28000/μ),Q(1)))
    require(all(x[1]<x[2] for x in results),"final closure inequality failed")
    return Dict(string(k)=>Dict("exact_value"=>string(v[1]),"exact_target"=>string(v[2])) for (k,v) in pairs(results))
end
