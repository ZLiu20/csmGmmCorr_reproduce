# Validation against Sun, McCaw and Lin's csmGmm package

The three scripts in this folder compare csmGmmCorr with the csmGmm package of Sun et al.
(https://github.com/ryanrsun/csmGmm). They are run by `make_all.sh` (skip them with `VALIDATION=0`).

| Script | Checks | Log |
|---|---|---|
| `validate_csmgmm.R` | csmGmm's EM started at our estimate (K = 4, Sigma = I; K = 2, fixed correlation); both fitted from the same start | `validate_csmgmm_log.txt` |
| `validate2_csmgmm.R` | A. one EM step from the same parameters; B. full fits from the same start on well-identified data (log-likelihood, lfdr, rejections); C. log-likelihood of both solutions on the sparse simulation data of `validate_csmgmm.R` | `validate2_log.txt` |
| `diag_incon.R` | our incongruity count against `csmGmm::check_incongruous`, and the diagnostics of Supplementary Section S4.6.3 | `diag_incon_log.txt` |
| `sim_convergence.R` | convergence of every csmGmm fit in the simulation studies (no replicate without a converged start) | `sim_convergence_by_design.csv` |

The numbers of Supplementary Sections 6.4 (comparison with csmGmm) and 7.3 come from
`validate2_log.txt` (parts A and B) and `diag_incon_log.txt`.

## Version of csmGmm

All scripts and `setup.R` install **csmGmm 0.5.0, GitHub commit `4a03af5`** (2026-06-26), so that
every run compares with the same code.

Changes from 0.4.0 (commit `a116d45`, 2025-09-16) to 0.5.0 in the functions used here
(`git diff a116d45 4a03af5`):

- `symm_fit_ind_EM`, `symm_fit_cor_EM`, `symm_fit_cor_EM_rho`: the default of `checkpoint` changed
  from `TRUE` to `FALSE`, and the progress output from `cat()` to `message()`. The EM iterations
  themselves are unchanged, so the fits are the same under both versions.
- `check_incongruous`: support for K = 4 was added (`find_4d`); in 0.4.0 the function stopped for
  K > 3. The K = 4 counts of Supplementary Section S4.6.3 therefore require 0.5.0 (or the upstream
  commit `c97c98c`, 2025-12-30, which already contained this change under the version number 0.4.0).

## Reading the first check of `validate_csmgmm.R`

On the sparse pleiotropy data of `validate_csmgmm.R` (few features with three or four nonzero
effects), csmGmm's EM started at our estimate can move far in the weights and magnitudes of the
nearly empty supports while the log-likelihood barely changes; part C of `validate2_csmgmm.R`
reports the log-likelihood of both solutions. Exact agreement is checked on well-identified data
in part B of `validate2_csmgmm.R`.
