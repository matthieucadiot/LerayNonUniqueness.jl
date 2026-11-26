# Only full, matching columns enter the final interval comparisons.
include("Setup.jl")
include("Files.jl")
par=residual_parameters();setprecision(par.prec)
data=load_data(joinpath(DATA,"Profile"))
nu=length(data.U)-1;nv=length(data.V)-1
signature=residual_signature()
S=zeros(Interval{BigFloat},19,nu)
for i in 1:nu
    S[:,i]=read_column(column_file(i),i,nv,par,signature).S
end
far=far_field_data(par.prec,data.λ̃)
residual=assemble(S[1:14,:],far,data,par)
verify_bounds(residual)
L=sqrt(U0(interval(BigFloat,data.L)))
ext=profile_extra_exterior(far,L)
extras=sum(S[15:17,:];dims=2)
gextras=(nRΨ2=extras[1],gΨ̃2_cells=extras[2],gΨ̃2_ext=ext.gΨ̃2_ext,RΩ̃2_cells=extras[3],RΩ̃2_ext=ext.RΩ̃2_ext)
c=profile_response_bounds(residual,gextras;complete=true)
exterior=table1_exterior(far,L)
table1=assemble_table1(S[18:19,:],c,exterior;complete=true)
verify_table1(table1)
residual_signature()==signature || error("the input or code changed during assembly")
file=joinpath(OUTPUT,"Residual_Table1.jld2")
mkpath(dirname(file));temporary,io=mktemp(dirname(file));close(io)
try
    jldsave(temporary;run_id=ENV["CLEAN_RUN_ID"],certified=true,residual=residual,table1=table1,profile_constants=c,profile_extra_integrals=gextras,
            exterior=exterior,S=S,columns=collect(1:nu),rows=collect(1:nv),par=par,signature=signature,
            profile_hash=signature["Profile/Interval_data.jld2"])
    Base.Filesystem.rename(temporary,file)
finally
    isfile(temporary)&&rm(temporary)
end
display((drift=table1.drift,test_fields=table1.test,profile_constants=c))
println("Every residual and Table 1 inequality is certified. Saved ",file)
