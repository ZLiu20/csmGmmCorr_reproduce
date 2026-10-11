// Genome-wide kernels for the real-data analysis (plain R C API; no Rcpp needed).
// Compile:  R CMD SHLIB -o rd_kernels.so rd_kernels.cpp      (done by setup.sh)
//
// Order of evidence (main text, Definition 1; csmGmmCorr::csm_adjust):
//   pleiotropy / mediation : class = open sign orthant of z (a zero coordinate -> own class),
//                            w = |z|;
//   replication            : class = directional branch, + when the (t+1)-th largest z is at
//                            least the (t+1)-th largest -z; w = z on +, -z on -.
// Feature j is at least as strong as feature i when they share the class and w_j >= w_i
// coordinatewise.
//
// rd_majorant   : monotone adjustment  v_i = max{ s_j : j at least as strong as i }  (exact,
//                 KD-tree with upper-orthant queries; equals csm_adjust(all = TRUE)).
// rd_incon_dec  : decision-level incongruity. For each non-rejected i, TRUE when some rejected
//                 j of the same class satisfies w_j <= w_i coordinatewise with w_j != w_i
//                 (i dominates a rejected feature but is not rejected).
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>

namespace rdk {

// class labels and transformed coordinates; z is column-major n x k
inline void represent(const double* z, int n, int k, int mode, int t,
                      std::vector<long long>& group, std::vector<double>& w) {
  group.assign(n, 0); w.assign((size_t)n * k, 0.0);
  std::vector<double> v(k);
  for (int i = 0; i < n; ++i) {
    if (mode == 0) {
      long long g = 0; bool zero = false;
      for (int d = 0; d < k; ++d) {
        double x = z[i + (size_t)d * n];
        if (x == 0.0) zero = true;
        if (x > 0.0) g |= (1LL << d);
        w[(size_t)i * k + d] = std::fabs(x);
      }
      group[i] = zero ? ((1LL << k) + i) : g;
    } else {
      int p = t + 1;
      for (int d = 0; d < k; ++d) v[d] = z[i + (size_t)d * n];
      std::sort(v.begin(), v.end());
      double plus = v[k - p], minus = -v[p - 1];
      double s = (plus >= minus) ? 1.0 : -1.0;
      group[i] = (plus >= minus) ? 0 : 1;
      for (int d = 0; d < k; ++d) w[(size_t)i * k + d] = s * z[i + (size_t)d * n];
    }
  }
}

struct Node { int b, e, l, r; double mx; std::vector<double> lo, hi; };

class Tree {
 public:
  Tree(const std::vector<double>& w, int k, const double* s, std::vector<int> rows, int leaf)
      : w_(w), k_(k), s_(s), rows_(rows), leaf_(std::max(2, leaf)) {
    nodes_.reserve(2 * (rows_.size() / leaf_ + 1));
    root_ = rows_.empty() ? -1 : build(0, (int)rows_.size());
  }
  // max score over rows with w_j >= q (all coordinates); returns best if none larger
  double qmax(const double* q, double best) const { return root_ < 0 ? best : qmax_(root_, q, best); }
  // exists row with w_j <= q (all) and w_j != q
  bool qdom(const double* q) const { return root_ >= 0 && qdom_(root_, q); }
 private:
  const std::vector<double>& w_; int k_; const double* s_; std::vector<int> rows_; int leaf_, root_;
  std::vector<Node> nodes_;
  inline double c(int row, int d) const { return w_[(size_t)row * k_ + d]; }
  int build(int b, int e) {
    int id = (int)nodes_.size();
    Node nd; nd.b = b; nd.e = e; nd.l = nd.r = -1; nd.mx = -std::numeric_limits<double>::infinity();
    nd.lo.assign(k_, std::numeric_limits<double>::infinity());
    nd.hi.assign(k_, -std::numeric_limits<double>::infinity());
    for (int p = b; p < e; ++p) {
      int row = rows_[p];
      if (s_) nd.mx = std::max(nd.mx, s_[row]);
      for (int d = 0; d < k_; ++d) { double x = c(row, d); nd.lo[d] = std::min(nd.lo[d], x); nd.hi[d] = std::max(nd.hi[d], x); }
    }
    nodes_.push_back(nd);
    if (e - b <= leaf_) return id;
    int sd = 0; double spread = -1;
    for (int d = 0; d < k_; ++d) { double sp = nodes_[id].hi[d] - nodes_[id].lo[d]; if (sp > spread) { spread = sp; sd = d; } }
    if (!(spread > 0)) return id;                       // identical rows: keep as a leaf
    int m = b + (e - b) / 2;
    std::nth_element(rows_.begin() + b, rows_.begin() + m, rows_.begin() + e,
                     [this, sd](int x, int y) { double a = c(x, sd), bb = c(y, sd); return a != bb ? a < bb : x < y; });
    int L = build(b, m), R = build(m, e);
    nodes_[id].l = L; nodes_[id].r = R;
    return id;
  }
  double qmax_(int id, const double* q, double best) const {
    const Node& nd = nodes_[id];
    if (nd.mx <= best) return best;
    bool all_in = true;
    for (int d = 0; d < k_; ++d) { if (nd.hi[d] < q[d]) return best; if (nd.lo[d] < q[d]) all_in = false; }
    if (all_in) return std::max(best, nd.mx);
    if (nd.l < 0) {
      for (int p = nd.b; p < nd.e; ++p) {
        int row = rows_[p]; if (s_[row] <= best) continue;
        bool ok = true; for (int d = 0; d < k_; ++d) if (c(row, d) < q[d]) { ok = false; break; }
        if (ok) best = s_[row];
      }
      return best;
    }
    int a = nd.l, b2 = nd.r; if (nodes_[b2].mx > nodes_[a].mx) std::swap(a, b2);
    best = qmax_(a, q, best);
    if (nodes_[b2].mx > best) best = qmax_(b2, q, best);
    return best;
  }
  bool qdom_(int id, const double* q) const {
    const Node& nd = nodes_[id];
    bool all_in = true, strict = false;
    for (int d = 0; d < k_; ++d) { if (nd.lo[d] > q[d]) return false; if (nd.hi[d] > q[d]) all_in = false; if (nd.hi[d] < q[d]) strict = true; }
    if (all_in && strict) return true;
    if (nd.l < 0) {
      for (int p = nd.b; p < nd.e; ++p) {
        int row = rows_[p]; bool ok = true, st = false;
        for (int d = 0; d < k_; ++d) { double x = c(row, d); if (x > q[d]) { ok = false; break; } if (x < q[d]) st = true; }
        if (ok && st) return true;
      }
      return false;
    }
    return qdom_(nd.l, q) || qdom_(nd.r, q);
  }
};

inline void majorant(const double* z, int n, int k, const double* s, int mode, int t, int leaf, double* out) {
  std::vector<long long> g; std::vector<double> w; represent(z, n, k, mode, t, g, w);
  std::map<long long, std::vector<int> > by; for (int i = 0; i < n; ++i) by[g[i]].push_back(i);
  for (int i = 0; i < n; ++i) out[i] = s[i];
  for (auto& it : by) {
    if (it.second.size() <= 1) continue;
    Tree tr(w, k, s, it.second, leaf);
    for (int i : it.second) out[i] = tr.qmax(&w[(size_t)i * k], s[i]);
  }
}

inline void incon_dec(const double* z, int n, int k, const int* rej, int mode, int t, int leaf, int* flag) {
  std::vector<long long> g; std::vector<double> w; represent(z, n, k, mode, t, g, w);
  std::map<long long, std::vector<int> > rj, nr;
  for (int i = 0; i < n; ++i) { flag[i] = 0; (rej[i] ? rj : nr)[g[i]].push_back(i); }
  for (auto& it : rj) {
    auto f = nr.find(it.first); if (f == nr.end()) continue;
    Tree tr(w, k, nullptr, it.second, leaf);
    for (int i : f->second) flag[i] = tr.qdom(&w[(size_t)i * k]) ? 1 : 0;
  }
}

}  // namespace rdk

#ifndef RD_TEST
#include <R.h>
#include <Rinternals.h>
extern "C" {
SEXP rd_majorant(SEXP z, SEXP s, SEXP mode, SEXP t, SEXP leaf) {
  int n = Rf_nrows(z), k = Rf_ncols(z);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  rdk::majorant(REAL(z), n, k, REAL(s), Rf_asInteger(mode), Rf_asInteger(t), Rf_asInteger(leaf), REAL(out));
  UNPROTECT(1); return out;
}
SEXP rd_incon_dec(SEXP z, SEXP rej, SEXP mode, SEXP t, SEXP leaf) {
  int n = Rf_nrows(z), k = Rf_ncols(z);
  SEXP out = PROTECT(Rf_allocVector(LGLSXP, n));
  rdk::incon_dec(REAL(z), n, k, LOGICAL(rej), Rf_asInteger(mode), Rf_asInteger(t), Rf_asInteger(leaf), LOGICAL(out));
  UNPROTECT(1); return out;
}
}
#endif
