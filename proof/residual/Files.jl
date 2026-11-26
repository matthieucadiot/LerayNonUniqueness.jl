# Completed columns and resumable row prefixes use the same nineteen quantities.
column_file(i)=joinpath(OUTPUT,"columns","column_$i.jld2")
partial_file(i)=joinpath(OUTPUT,"columns","partial_$i.jld2")

function valid_quantities(S,prec)
    size(S)==(19,) || return false
    return all(x->x isa Interval{BigFloat}&&isguaranteed(x)&&isfinite(inf(x))&&isfinite(sup(x))&&
                  inf(x)≥0&&precision(inf(x))==precision(sup(x))==prec,S)
end

function save_column(file,S,i,nv,done,par,signature)
    i≥1&&nv≥1&&0≤done≤nv || error("invalid column coverage")
    valid_quantities(S,par.prec) || error("invalid column quantities")
    mkpath(dirname(file));temporary,io=mktemp(dirname(file));close(io)
    try
        jldsave(temporary;run_id=ENV["CLEAN_RUN_ID"],schema=1,S=S,column=i,nv=nv,rows_done=done,rows=collect(1:done),par=par,signature=signature)
        Base.Filesystem.rename(temporary,file)
    finally
        isfile(temporary)&&rm(temporary)
    end
end

function read_column(file,i,nv,par,signature;complete=true)
    r=load(file)
    r["run_id"]==ENV["CLEAN_RUN_ID"] || error("Residual column belongs to a different run: $file")
    get(r,"schema",0)==1&&r["column"]==i&&r["nv"]==nv&&r["par"]==par&&r["signature"]==signature || error("$file does not match this run")
    done=r["rows_done"]
    done isa Integer&&0≤done≤nv&&r["rows"]==collect(1:done) || error("$file has invalid row coverage")
    (complete ? done==nv : done<nv) || error("$file has incomplete or unexpected coverage")
    valid_quantities(r["S"],par.prec) || error("$file has invalid interval quantities")
    return (S=r["S"],rows_done=Int(done))
end

function accumulate!(S,b)
    length(S)==length(b)==19 || error("a profile cell has nineteen quantities")
    S[1:10]=S[1:10]+b[1:10]
    S[11:14]=max.(S[11:14],b[11:14])
    S[15:18]=S[15:18]+b[15:18]
    S[19]=max(S[19],b[19])
    return S
end
