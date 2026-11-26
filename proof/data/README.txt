FIXED INTERVAL INPUT DATA

Profile/Interval_data.jld2
  Four profile spline coefficient arrays, knots, original partitions, cutoff
  parameters, and the interval containing the exact decimal eigenvalue.

responses_odd.jld2, responses_even.jld2
  The fixed response spline numerators, knots and weights. The even file also
  contains its scalar coefficients and four exterior harmonic coefficients.

functionals_odd.jld2, functionals_even.jld2
  The selected functional directions, coefficients, retained modes, bumps and
  supplied normalized eigenfunction approximations.

Proposition_4_21_inputs.jld2
  The full cell construction inputs, including all nine spectral bases, all
  207 nonconstant modes and the selected directions. Boxes and the remaining
  columns of the matrices X are constructed during the proof.

Eigenfunctions/*.jld2
  The 207 interval eigenpair approximations for the recurrence and matching
  proof of Appendix A.

Eigenfunctions_normalized.jld2
  The fixed normalized coefficients whose errors Appendix A must certify.
  They agree with the approximations used by the functionals and by the cell
  construction. This file supplies coefficients, not an assumed proof result.

All proof inputs are already JLD2 interval data; no Python conversion or
floating-point eigensolver is needed for a run. Integers, names and dimensions
are metadata. The workers validate interval precision, finiteness, parameters,
coverage and source/input hashes. Input data are kept separate from results.

SHA256SUMS records the fixed inputs, including the weighted response data. Each code group also has a manifest
for the subset it needs. Identical shared files in different Della archives
have identical checksums.

FLOAT64 COERCIVITY DATA
Coercivity/Profile_float64.jld2 and Coercivity/Proposition_float64.jld2 store
outward Float64 interval coefficients for the new coercivity version.
Small exact topology and stable seed parameters retain their source precision.
The conversion checked every original coefficient and records original-input
and conversion-code fingerprints. Original residual and eigenfunction files
are retained unchanged.

CLEAN SELECTIVE-WEIGHT INPUTS

Coercivity/selected_factors.tsv gives the exact rational choices used in the
paper. Coercivity/Fixed_weighted_responses.jld2 embeds the corresponding odd
spline numerators, all 58 adjusted matrix penalties and the corrected even
scalar coefficients. The original responses and selected functional sources
above remain available to audit the construction. The table2/clean_prepare.jl
builder and table2/clean_inputs.jl reader are under ../coercivity. Final results
from a new run go to its separate scratch directory.
