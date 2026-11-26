# LerayNonUniqueness.jl

Code, spline data and certification records for *On the non-uniqueness of solutions of the axi-symmetric swirl-free Navier–Stokes equations, II*, by Matthieu Cadiot, Alexandru Ionescu, Hao Jia and Stan Palasek. [Paper on arXiv](https://arxiv.org/abs/2611.27984) · [Versioned downloads](https://github.com/matthieucadiot/LerayNonUniqueness.jl/releases/tag/paper-v1-2611)

The clean computation finished on 26 November 2026 with Julia 1.12.3. It certified the numerical bounds used in the paper, including all 58 responses, the original Gershgorin inequalities for all 15 matrix cases and the final contraction estimates. The response sweep covered 57,725,138 parent cells in 67,682 columns. The final verifier reports `all_checks_passed = true`. [Results for Tables 1–5](TABLES.md) · [Original final report](verification.toml)

## Choose a download

| What you want to do | Download | Compressed size |
| --- | --- | ---: |
| Check the recorded results for every table | [Certification evidence](https://github.com/matthieucadiot/LerayNonUniqueness.jl/releases/download/paper-v1-2611/certification-evidence.zip) | 73.9 MiB |
| Obtain the fixed spline and spectral inputs or rerun the complete proof | [Clean proof code and inputs](https://github.com/matthieucadiot/LerayNonUniqueness.jl/releases/download/paper-v1-2611/Proof_Float64_Clean_2611.zip) | 321.2 MiB |
| Inspect the four approximate profile fields in NumPy or reproduce the nonlinear and eigenpair solves | [Profile construction](https://github.com/matthieucadiot/LerayNonUniqueness.jl/releases/download/paper-v1-2611/profile-construction.zip) | 41.6 MiB |
| Reassemble the earlier 6 November residual calculation from its 804 original column records | [Historical residual certification](https://github.com/matthieucadiot/LerayNonUniqueness.jl/releases/download/paper-v1-2611/certification-1.zip) | 4.5 MiB |

The evidence download stands alone for the Python check. The clean proof download stands alone for recomputation once the pinned software is installed. For the optional historical replay, also extract the clean proof download so that both calculations can use the same interval spline file. [Archive checksums](SHA256SUMS)

## Check the evidence

With Python 3.11 or later, run the following from the directory containing the downloaded evidence ZIP.

```sh
unzip -q certification-evidence.zip
python3 certification-evidence/verify_evidence.py
```

The checker verifies the final report against its recorded SHA-256, checks the original result files linked from that report, and compares the recorded response bounds, matrix margins and rational closing estimates with their targets. The Julia table reader exposes the interval values in the underlying certificates. To inspect the numerical methods in a browser, start with [the final acceptance code](proof/scripts/clean_final.jl). For the Julia environments, a small numerical example, and instructions for a fresh full run, see [the quickstart](QUICKSTART.md).

## Locate the splines

Within the clean proof ZIP, `data/Profile/Interval_data.jld2` contains the knots and the four coefficient arrays for the approximate stationary profile and its unstable eigenmode. The optional profile construction ZIP provides the original double-double coefficients in two NumPy files, together with the numerical solver, seed and interval-conversion code. Rebuilding the interval file from those NumPy outputs with Julia 1.12.3 reproduced every numerical byte of the published input. The only byte that differed in the file was the Julia patch version in the writer’s header (1.12.4 in the original).

The auxiliary Table 2 response splines are in `data/responses_odd.jld2` and `data/responses_even.jld2`. The proof uses `data/Coercivity/Fixed_weighted_responses.jld2` together with the exact selective factors and the fixed input data. The clean code checks the relationship between the supplied responses and their weighted form. The profile-construction solver generates the stationary profile and its unstable eigenpair. The 58 auxiliary responses are supplied as fixed inputs for verification. [Input inventory](proof/data/README.txt) · [Selective factors](proof/data/Coercivity/selected_factors.tsv)

The final run identity is `2d9cec66690641fa88a3aac70215c5cb`. Its verification report has SHA-256 `42fc3d336755a2c10b9d1338ba298e66bc86b7b1cbe8a3be32451b2703292104`.
