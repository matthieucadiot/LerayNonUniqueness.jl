# Profile_response_bounds.jl : the profile constants used in Lemma 6.2 and the even response projection.
# These are derived from computed integrals, not copied from the certification targets of Table 1.

profile_response_nonnegative(x) = x isa Interval && isguaranteed(x) && isfinite(inf(x)) && isfinite(sup(x)) && inf(x)≥0

# b is the result of the residual assemble(). In addition to its 14 cell quantities, the shared profile
# sweep supplies ||R_Psi||_2², ||grad Psi_tilde5||_2² and ||R Omega_tilde||_2².
# The last two have separate cell and spherical-exterior contributions. The exterior contributions are
# upper majorants only; their lower endpoints must not enter a lower bound on the actual norm.
function profile_response_bounds(b,gextras;complete=false)
    values=(b.wRΨin,b.wRΨband,b.rRΨ̃,b.nΩ̃,b.sΩ,b.sΩfar,b.IoΩ,b.ℛ̃,
            gextras.nRΨ2,gextras.gΨ̃2_cells,gextras.gΨ̃2_ext,gextras.RΩ̃2_cells,gextras.RΩ̃2_ext)
    all(profile_response_nonnegative,values) || error("profile response constants require nonnegative guaranteed finite integrals and bounds")
    RΨZ=sqrt(gextras.nRΨ2+b.wRΨin*b.wRΨin+b.wRΨband*b.wRΨband)
    gΨ̃5=sqrt(gextras.gΨ̃2_cells+I"[0,1]"*gextras.gΨ̃2_ext)
    lower=sqrt(gextras.gΨ̃2_cells)-b.rRΨ̃
    inf(lower)>0 || error("the computed corrected even potential has no positive gradient-norm lower bound")
    gΨ̃_min=interval(Float64,inf(lower))                               # an exactly represented downward-rounded lower bound
    gΨ̃=gΨ̃5+b.rRΨ̃                                                    # upper majorant for ||grad Delta^-1 Omega_tilde||_2
    nΩ̃Z=sqrt(b.nΩ̃*b.nΩ̃+gΨ̃*gΨ̃)
    RΩ̃=sqrt(gextras.RΩ̃2_cells+I"[0,1]"*gextras.RΩ̃2_ext)
    sSΩ̃=RΩ̃+I"5"*gΨ̃
    all(profile_response_nonnegative,(RΨZ,gΨ̃_min,gΨ̃,nΩ̃Z,RΩ̃,sSΩ̃)) || error("invalid assembled profile constants")
    return (Cb=I"8"/I"100",RΨZ=RΨZ,sΩ=b.sΩ,IoΩ=b.IoΩ,sΩfar=b.sΩfar,rRΨ̃=b.rRΨ̃,
            gΨ̃_min=gΨ̃_min,nΩ̃Z=nΩ̃Z,sSΩ̃=sSΩ̃,ℛ̃=b.ℛ̃,δ=I"1"/I"1e6",
            gΨ̃5=gΨ̃5,gΨ̃=gΨ̃,RΩ̃=RΩ̃,certified=complete===true)
end

function assemble_profile_response(S14,gextras,far,data,par;complete=false)
    b=assemble(S14,far,data,par)
    return profile_response_bounds(b,gextras;complete=complete)
end
