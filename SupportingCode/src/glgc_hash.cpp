#include <cstdio>
#include <string>
#include "blake2b8.h"
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

extern "C" SEXP glgc_blake2b8_hex(SEXP keys) {
  const R_xlen_t n = XLENGTH(keys);
  SEXP answer = PROTECT(Rf_allocVector(STRSXP, n));
  char buffer[17];
  for (R_xlen_t i = 0; i < n; ++i) {
    const char *value = Rf_translateCharUTF8(STRING_ELT(keys, i));
    const unsigned long long digest = static_cast<unsigned long long>(blake2b8(std::string(value)));
    std::snprintf(buffer, sizeof(buffer), "%016llx", digest);
    SET_STRING_ELT(answer, i, Rf_mkChar(buffer));
  }
  UNPROTECT(1);
  return answer;
}
