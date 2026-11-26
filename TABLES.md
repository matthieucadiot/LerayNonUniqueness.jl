# Certification evidence for Tables 1–5

The clean verification finished on Della on **26 November 2026 at 08:11 UTC**. The final receipt records `all_checks_passed = true`. The original job number is `18207736`. The SHA-256 of `verification.toml` is `42fc3d336755a2c10b9d1338ba298e66bc86b7b1cbe8a3be32451b2703292104`.

The table numbers and page references below follow the paper on arXiv (v1). The numerical columns in this guide come from the original aggregate certificates. For display, upper bounds are rounded up and lower bounds down. The reader compares the full stored endpoints with exact decimal targets. [Paper, arXiv:2611.27984](https://arxiv.org/abs/2611.27984v1)

| Paper | What was checked | Original records inside the evidence ZIP |
| --- | --- | --- |
| Table 1, p. 10 | Bounds for the nonlinear solution and the eigenpair, including both fixed-point residuals | `results/residual/Residual_Table1.jld2` |
| Table 2, p. 21 | Four sums for each parity, using 26 odd and 32 even functionals with the weights used in the paper | `results/coercivity/Table2_bounds.jld2`, `results/coercivity/table2_checks.toml` |
| Table 3, p. 35 | The 13 constituent residual and coefficient bounds, including the exterior region | `results/residual/Residual_Table1.jld2`, saved field `residual` |
| Table 4, p. 52 | Cell parameters, numbers of selected directions, complementary gaps and matrix inequalities | `results/coercivity/proposition/` contains `Bounds.jld2`, `Boxes.jld2`, `table4_checks.toml` and 15 files in `matrices/` |
| Table 5, p. 78 | Weight parameters and highest certified mode for each of nine weights, with an enclosure for every nonconstant eigenpair | `results/eigenfunctions/Eigenfunctions.jld2`, saved fields `certificates` and `constants` |

## Read and verify the results

The numerical summaries below can be read directly on GitHub. For a local check, extract the evidence and clean proof ZIPs together. Starting in the directory that contains the two extracted folders, use Python 3.11 or later and Julia 1.12.3. The first Julia command installs the dependencies pinned by the original proof.

```sh
EVIDENCE="$PWD/certification-evidence"
PROOF="$PWD/Proof_Float64_Clean_2611"
python3 "$EVIDENCE/verify_evidence.py" "$EVIDENCE" --proof "$PROOF"
julia --startup-file=no --project="$PROOF/residual" -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project="$PROOF/residual" "$EVIDENCE/read_tables.jl" "$EVIDENCE" "$PROOF"
```

Python checks the receipt, digests, reported bounds and exact closing comparisons. Julia reads the interval objects, checks the paper's numerical inequalities, checks the row definitions against the sealed proof input and prints the five tables. A full independent recomputation starts from the fixed spline and coefficient files in the clean proof ZIP. Its entry points and Slurm instructions are in that ZIP's README.

## Table 1. Quantitative bounds for Proposition 2.1

The two residual entries are the fully assembled norms, including exterior contributions. The odd upper bound also satisfies the sharper 5.85 × 10⁻¹² estimate at the end of Section 3.3. The aggregate records all 804 × 575 profile cells.

| Quantity | Certified bound | Paper |
| --- | --- | --- |
| Nonlinear residual in Z′ | ≤ 5.65272913e-12 | ≤ 6e-12 |
| Eigenpair residual in Z′ | ≤ 1.35562331e-12 | ≤ 1.5e-12 |
| Stream residual in Z | ≤ 1.18054384e-9 | ≤ 1.251e-9 |
| Even potential gradient in L² | [0.883628755, 0.883628756] | [0.88, 0.89] |
| rΩ* in L∞ | ≤ 12.4765492 | ≤ 13 |
| Weighted even profile in L² | ≤ 22.4987168 | ≤ 22.55 |
| Drift field in L∞ | ≤ 22.2970670 | ≤ 22.45 |
| Joint L² norm of the two test fields | ≤ 14.3261332 | ≤ 14.4 |

## Table 2. Functional and response sums

Each entry is a sum over one parity class. The response norm in the last product is Z for odd responses and Z₀ for even responses. The run covered 67,682 columns, with 57,725,138 parent cells and 57,901,397 evaluated cells under the subdivision rule. That rule selected 4,289 columns. The two residual sums and two response products also meet the sharper printed inequalities in (6.65), respectively 4.790 × 10⁻⁴, 7.131 × 10⁻⁵, 15.244 and 141.786.

| Parity and sum | Certified upper | Paper upper |
| --- | --- | --- |
| odd functional norms | 18.8767713 | 19.2 |
| odd squared functional norms | 19.5851833 | 19.86 |
| odd functional × residual norms | 4.78997388e-4 | 5e-4 |
| odd response × functional norms | 15.2437228 | 15.25 |
| even functional norms | 33.7913575 | 34 |
| even squared functional norms | 60.2070662 | 61 |
| even functional × residual norms | 7.13001591e-5 | 1e-3 |
| even response × functional norms | 141.785109 | 146 |

## Table 3. Bounds on residuals and coefficients

The quantities are listed separately for each entry printed on a single manuscript row. The saved Table 3 fields and the two assembled residuals are checked in `residual/Norms.jl`. The other profile bounds in Table 1 are checked in `residual/Table1_fields.jl`.

| Quantity | Certified bound | Paper |
| --- | --- | --- |
| rΩ* in L∞ | ≤ 12.4765492 | ≤ 13 |
| rΩ̃* in L∞ | ≤ 0.0240710803 | ≤ 0.025 |
| rΩ* in L∞ outside B₉₀₀ | ≤ 4.97751905e-4 | ≤ 5e-4 |
| rΩ̃* in L∞ outside B₉₀₀ | ≤ 7.80871864e-19 | ≤ 1e-18 |
| Iₒ(Ω*, 900) | ≤ 0.370775068 | ≤ 0.38 |
| Iₒ(Ω̃*, 900) | ≤ 3.03105715e-5 | ≤ 3.2e-5 |
| Ω̃* in L² | [1.00499999, 1.00500000] | [1, 1.01] |
| √w times the Ω residual in L² | ≤ 4.20802982e-12 | ≤ 4.3e-12 |
| r times the Ω̃ residual in L² | ≤ 1.16563655e-12 | ≤ 1.25e-12 |
| √w times the interior Ψ residual in L² | ≤ 6.60939856e-14 | ≤ 6.8e-14 |
| r times the Ψ̃ residual in L² | ≤ 1.50986956e-14 | ≤ 1.6e-14 |
| Band Ψ residual in L² | ≤ 5.33972597e-13 | ≤ 5.6e-13 |
| √w times the band Ψ residual in L² | ≤ 1.18054371e-9 | ≤ 1.25e-9 |

## Table 4. Cell construction and matrix certificates

Wᵢ, Λ* and nᵢ are the printed inputs. The last three columns are interval bounds obtained for that row. The gap is the lower bound for the omitted modes. The last column is the smallest margin over every selected and unselected direction with the prescribed weights. The 15 matrices pass the original Gershgorin inequalities, the norm threshold 30, the orthogonality-defect threshold 10⁻¹² and the loss threshold 9 × 10⁻¹¹. The inputs select 26 odd vector directions, 23 even vector directions and 9 even scalar directions. The final verifier also checked the 2,482 box columns, 804 × 575 far-field cells, 2,521 × 321 support cells and the analytic exterior bounds.

| Row | Wᵢ | Λ* | nᵢ | Verified gap ≥ | Norm ≤ | Least Gershgorin margin ≥ |
| --- | --- | --- | --- | --- | --- | --- |
| odd_T0U0 | 3.76 | 12.5455 | 9 | 12.5455500224034 | 11.8911420 | 2.17892714e-6 |
| odd_T0U1 | 0.504 | 1.7237 | 1 | 1.72372778963434 | 1.74830053 | 6.57364238e-4 |
| odd_T1U0 | 4.96 | 16.6307 | 11 | 16.6307636012840 | 15.7292970 | 2.58947836e-5 |
| odd_T1U1 | 0.489 | 1.6562 | 3 | 1.65627166177675 | 1.73589143 | 4.71938954e-6 |
| odd_T2U0 | 0.362 | 1.3053 | 1 | 1.30532846537771 | 1.16216112 | 5.61949779e-5 |
| odd_T2U1 | 0.235 | 0.854 | 1 | 0.854032359201908 | 0.818864610 | 2.97064937e-5 |
| even_T0U0 | 3.71 | 12.3774608652 | 12 | 12.3774608652476 | 11.7382762 | 9.26838008e-4 |
| even_T0U1 | 0.494 | 1.7237277896 | 0 | 1.72372778963434 | 1.76125820 | 0.00272353774 |
| even_T1U0 | 4.72 | 15.7151282883 | 9 | 15.7151282883570 | 15.0747557 | 7.47770326e-4 |
| even_T1U1 | 0.477 | 1.623936892 | 2 | 1.62393689203068 | 1.65788736 | 9.99989230e-4 |
| even_T2U0 | 0.285 | 1.0866964506 | 0 | 1.08669645068305 | 0.990485480 | 0.0155007452 |
| even_T2U1 | 0.184 | 0.6501580511 | 0 | 0.650158051176523 | 0.728558074 | 0.0113289904 |
| scalar_S0V0 | 0.706 | 11.2020998964 | 6 | 11.2020998964803 | 2.12397352 | 9.38927927e-4 |
| scalar_S1V0 | 0.706 | 11.2110651424 | 3 | 11.2110651424916 | 2.11210500 | 2.77541368e-4 |
| scalar_S2V1 | 0.00238 | 0.0378391209 | 0 | 0.0378391209721839 | 0.00690497432 | 2.26407923e-6 |

## Table 5. Spectral weights and verified modes

The radial weights T₀, T₁, T₂, S₀ and S₂ carry the factor r³. T₁ and S₁ give the same weight. For each of the nine weights the zero eigenvalue has an exact constant mode. The aggregate encloses each nonconstant mode from 1 through the listed highest index and records its zero count. The degree is 200 and the working precision is 128 bits. Across the 207 nonconstant modes, the eigenvalue error is at most 3.3 × 10⁻³¹, the uniform error in derivative orders 0, 1 and 2 is at most 1.1 × 10⁻²⁴, and the contraction bound is at most 4.1 × 10⁻¹².

| Weight | (a, x₂, x₃, b) | Highest verified index | Largest eigenvalue error ≤ | Largest error in orders 0–2 ≤ |
| --- | --- | --- | --- | --- |
| T0 | (0.5, 1, 15, 17) | 19 | 1.06052554e-33 | 9.26187622e-29 |
| T1 | (15, 17, 29, 31) | 21 | 2.01204311e-31 | 6.89112042e-29 |
| T2 | (29, 31, 42, 45) | 6 | 3.28639155e-31 | 9.60827607e-30 |
| S0 | (0.1, 0.2, 15, 17) | 18 | 5.47113500e-34 | 5.58153913e-27 |
| S2 | (29, 31, 152, 160) | 9 | 1.25662908e-32 | 3.70387922e-31 |
| U0 | (-7, -5, 5, 7) | 18 | 6.03551273e-33 | 8.23763192e-29 |
| U1 | (5, 7, 20, 23) | 8 | 1.24013994e-33 | 3.80319296e-31 |
| V0 | (-40, -35, 35, 40) | 84 | 1.25895054e-31 | 1.04919510e-24 |
| V1 | (-200, -180, 180, 200) | 24 | 8.72267937e-34 | 3.37990360e-32 |

## Closing estimates

The seven exact rational comparisons are in the `closure` section of `verification.toml`. They cover the odd inverse bounds in Z and Z₀ (both below 164), the images of the odd and even fixed-point balls and both contraction constants. The final clean proof entry point is `scripts/clean_final.jl`. The original Della execution is recorded in `logs/final.01_18207736_4294967294.out`.
