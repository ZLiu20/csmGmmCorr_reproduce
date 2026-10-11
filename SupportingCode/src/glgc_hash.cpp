// BLAKE2b-64 hash of variant keys, for the deterministic choice of one fitting variant per
// 10-kb window in Data/prepare_glgc.R (identical to Python hashlib.blake2b(key, digest_size=8)).
#include <cstdio>
#include "blake2b8.h"
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

extern "C" {
// Returns the 64-bit digests as 16-character hexadecimal strings (big-endian), so that the
// order of the strings equals the order of the unsigned integers.
SEXP glgc_blake2b8_hex(SEXP keys) {
  R_xlen_t n = XLENGTH(keys);
  SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
  char buf[17];
  for (R_xlen_t i = 0; i < n; i++) {
    const char* s = Rf_translateCharUTF8(STRING_ELT(keys, i));
    unsigned long long h = (unsigned long long)blake2b8(std::string(s));
    std::snprintf(buf, sizeof(buf), "%016llx", h);
    SET_STRING_ELT(out, i, Rf_mkChar(buf));
  }
  UNPROTECT(1);
  return out;
}
}
