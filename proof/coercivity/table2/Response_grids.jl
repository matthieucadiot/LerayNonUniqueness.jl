# Pure grid selection shared by response workers and work-list construction.
# response_partition is supplied by the including namespace. Keep all supplied
# spectral junctions in both grids; only the opposite response knots are omitted.
function response_grids(profile,odd,even,functionals=())
    return (odd=response_partition(profile,odd,odd,functionals),
            even=response_partition(profile,even,even,functionals))
end

response_profile(t,odd::Bool) = odd ? t.odd_profile : t.even_profile
