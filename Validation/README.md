# Validation against csmGmm

The scripts in this folder compare correlated csmGmm calculations with the
independent csmGmm implementation of Sun, McCaw and Lin whenever that package
is installed. `validate_csmgmm.R` checks fitted local false-discovery rates from
the same starting point; `validate2_csmgmm.R` checks a single update and a
well-identified full fit; `diag_incon.R` compares ordering diagnostics.

The cross-package validation is deliberately conditional on the external
package version. Absence of csmGmm is reported as a skipped comparison, not as
evidence of equivalence. The package's own likelihood, posterior, constraint,
thread-equivalence, and decision identities are checked separately by
`csmGmmCorr/tests/run-tests.R`.

Set `VAL_N`, `VAL_TOL`, and `VAL_MAPS` to control the comparison size.  The
ordering study uses `SIM_N`, `DIAG_MAX_ITER`, `INCON_MAX`, and `CORES`.  The
small defaults exported by `tests/test_all.sh local` are smoke settings; paper
scale diagnostics must be run explicitly on suitable compute resources.
