# Response_assemble.jl : Lemma 6.2, Remark 6.3, the even projection, and the four rows of Table 2.
# All numerical entries are enclosures of upper majorants: only their upper endpoints assert norm bounds.
# A positive lower endpoint for such a majorant is not a lower bound for the actual norm.

response_nonnegative(x) = x isa Interval && isguaranteed(x) && isfinite(inf(x)) && isfinite(sup(x)) && inf(x)≥0

function check_response_inputs(xs)
    all(response_nonnegative,xs) || error("response assembly requires nonnegative guaranteed finite intervals")
end

# b = (nΦ,gΨ,wE,wD_in,D_band,wD_band,complete): norms of the raw odd response, its potential,
# the local residual, and the interior/band potential defects. qerror bounds ||r(q-qnum)||_2.
# c contains explicit profile bounds; its certified flag records whether these hypotheses were independently proved.
function response_bound_odd(b,qnum,qerror,c; functional_certified=false)
    check_response_inputs((b.nΦ,b.gΨ,b.wE,b.wD_in,b.D_band,b.wD_band,qnum,qerror,c.Cb,c.RΨZ,c.sΩ,c.IoΩ,c.sΩfar))
    d=sqrt(b.wD_in*b.wD_in+b.wD_band*b.wD_band)
    g=b.gΨ+d
    Φ=sqrt(b.nΦ*b.nΦ+g*g)
    coupling=c.sΩ*b.wD_in+sqrt(c.IoΩ*c.IoΩ*b.D_band*b.D_band+c.sΩfar*c.sΩfar*b.wD_band*b.wD_band)
    e=b.wE+coupling+c.Cb*c.RΨZ*b.nΦ+qerror
    q=qnum+qerror
    check_response_inputs((g,Φ,e,q))
    certified=get(b,:complete,false)===true && get(c,:certified,false)===true && functional_certified===true
    return (q=q,e=e,Φ=Φ,ΦZ=Φ,g=g,m0=zero(g),certified=certified)
end

# b = (nΦ,RΦ,gΨ,rE,rD,orth,complete). RΦ means ||sqrt(r²+z²) Phi||_2, not ||r Phi||_2.
# orth encloses the signed integral <Phi,Psi_tilde5>; it is measured, rather than assumed to vanish exactly.
# The output Phi is the norm of P0 Phi in Z0, including its projection correction.
function response_bound_even(b,qnum,qerror,c; functional_certified=false)
    orth=abs(b.orth)
    check_response_inputs((b.nΦ,b.RΦ,b.gΨ,b.rE,b.rD,orth,qnum,qerror,c.Cb,c.RΨZ,c.sΩ,c.rRΨ̃,c.gΨ̃_min,c.nΩ̃Z,c.sSΩ̃,c.ℛ̃,c.δ))
    inf(c.gΨ̃_min)>0 || error("the corrected even profile potential needs a positive gradient-norm lower bound")
    g=b.gΨ+b.rD
    ΦZ=sqrt(b.nΦ*b.nΦ+g*g)
    SΦ=b.RΦ+I"5"*g
    m0=(orth+g*c.rRΨ̃)/(c.gΨ̃_min*c.gΨ̃_min)
    Φ=ΦZ+c.δ*SΦ+m0*(c.nΩ̃Z+c.δ*c.sSΩ̃)
    e=b.rE+c.sΩ*b.rD+c.Cb*c.RΨZ*b.nΦ+qerror+m0*c.ℛ̃
    q=qnum+qerror
    check_response_inputs((g,ΦZ,Φ,m0,e,q))
    certified=get(b,:complete,false)===true && get(c,:certified,false)===true && functional_certified===true
    return (q=q,e=e,Φ=Φ,ΦZ=ΦZ,g=g,m0=m0,certified=certified)
end

function response_sums(B)
    isempty(B) && error("there are no responses to assemble")
    return (q=sum(b.q for b in B),q2=sum(b.q*b.q for b in B),
            qe=sum(b.q*b.e for b in B),Φq=sum(b.Φ*b.q for b in B))
end

# Q entries are (norm,error,certified). The error is mandatory: it may be zero only when exact equality
# or an equivalent direct enclosure of the exact functional has been established.
function assemble_responses(Bodd,Beven,Qodd,Qeven,c)
    length(Bodd)==length(Qodd) && length(Beven)==length(Qeven) || error("functional and response counts differ")
    odd=[response_bound_odd(Bodd[j],Qodd[j].norm,Qodd[j].error,c;functional_certified=Qodd[j].certified) for j in eachindex(Bodd)]
    even=[response_bound_even(Beven[j],Qeven[j].norm,Qeven[j].error,c;functional_certified=Qeven[j].certified) for j in eachindex(Beven)]
    return (odd=response_sums(odd),even=response_sums(even),odd_responses=odd,even_responses=even,
            n_odd=length(odd),n_even=length(even),δ=c.δ,certified=all(b.certified for b in odd)&&all(b.certified for b in even))
end

function table2_targets()
    return (odd=(q=I"192"/I"10",q2=I"1986"/I"100",qe=I"1"/I"2000",Φq=I"61"/I"4"),
            even=(q=I"34",q2=I"61",qe=I"1"/I"1000",Φq=I"146"))
end

# This is a threshold diagnostic even while external certificates or full-grid coverage are pending.
function table2_checks(t)
    targets=table2_targets()
    passed(x,y)=response_nonnegative(x) && sup(x)≤inf(y)
    odd=(q=passed(t.odd.q,targets.odd.q),q2=passed(t.odd.q2,targets.odd.q2),
         qe=passed(t.odd.qe,targets.odd.qe),Φq=passed(t.odd.Φq,targets.odd.Φq))
    even=(q=passed(t.even.q,targets.even.q),q2=passed(t.even.q2,targets.even.q2),
          qe=passed(t.even.qe,targets.even.qe),Φq=passed(t.even.Φq,targets.even.Φq))
    return (odd=odd,even=even,counts=t.n_odd==26&&t.n_even==32)
end

function verify_table2(t)
    t.certified || error("Table 2 is conditional: full coverage or a profile/functional certificate is missing")
    in_interval(big(1)//big(10)^6,t.δ) || error("Table 2 uses the exact value delta=1e-6 in the even Z0 norm")
    checks=table2_checks(t)
    checks.counts || error("Table 2 requires exactly 26 odd and 32 even responses")
    all(values(checks.odd))&&all(values(checks.even)) || error("a Table 2 target has not been verified: $checks")
    return true
end
