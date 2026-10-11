#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

namespace {

struct Representation {
  std::vector<long long> group;
  std::vector<double> coordinate;
};

int mode_value(SEXP mode) {
  if (TYPEOF(mode) == STRSXP) {
    std::string x = CHAR(STRING_ELT(mode, 0));
    return x == "replication" ? 1 : 0;
  }
  return Rf_asInteger(mode);
}

Representation represent(const double *z, int n, int k, int mode, int t) {
  Representation out;
  out.group.resize(n);
  out.coordinate.resize(static_cast<size_t>(n) * k);
  if (mode == 1) {
    const int rank = t;
    std::vector<double> sorted(k);
    for (int i = 0; i < n; ++i) {
      for (int d = 0; d < k; ++d) sorted[d] = z[i + static_cast<size_t>(d) * n];
      std::sort(sorted.begin(), sorted.end());
      const double lower = sorted[rank];
      const double upper = sorted[k - rank - 1];
      const double sign = upper >= -lower ? 1.0 : -1.0;
      out.group[i] = sign > 0 ? 1 : 0;
      for (int d = 0; d < k; ++d) out.coordinate[static_cast<size_t>(i) * k + d] = sign * z[i + static_cast<size_t>(d) * n];
    }
  } else {
    for (int i = 0; i < n; ++i) {
      long long group = 0;
      bool zero = false;
      for (int d = 0; d < k; ++d) {
        const double x = z[i + static_cast<size_t>(d) * n];
        if (x == 0) zero = true;
        if (x > 0) group |= (1LL << d);
        out.coordinate[static_cast<size_t>(i) * k + d] = std::fabs(x);
      }
      out.group[i] = zero ? ((1LL << k) + i) : group;
    }
  }
  return out;
}

struct Node {
  int begin, end, left, right;
  double maximum;
  std::vector<double> low, high;
};

class DominanceTree {
 public:
  DominanceTree(const std::vector<double> &x, int dimension,
                const double *score, const std::vector<int> &members, int leaf_size)
    : x_(x), dimension_(dimension), score_(score), order_(members), leaf_size_(std::max(2, leaf_size)) {
    root_ = order_.empty() ? -1 : build(0, static_cast<int>(order_.size()));
  }

  double upper_max(const double *query, double initial) const {
    return root_ < 0 ? initial : upper_max_node(root_, query, initial);
  }

  bool lower_strict_exists(const double *query) const {
    return root_ >= 0 && lower_exists_node(root_, query);
  }

 private:
  const std::vector<double> &x_;
  int dimension_;
  const double *score_;
  std::vector<int> order_;
  int leaf_size_;
  int root_;
  std::vector<Node> nodes_;

  double at(int row, int d) const { return x_[static_cast<size_t>(row) * dimension_ + d]; }

  int build(int begin, int end) {
    Node node;
    node.begin = begin; node.end = end; node.left = node.right = -1;
    node.maximum = -std::numeric_limits<double>::infinity();
    node.low.assign(dimension_, std::numeric_limits<double>::infinity());
    node.high.assign(dimension_, -std::numeric_limits<double>::infinity());
    for (int p = begin; p < end; ++p) {
      int row = order_[p];
      if (score_) node.maximum = std::max(node.maximum, score_[row]);
      for (int d = 0; d < dimension_; ++d) {
        double value = at(row, d);
        node.low[d] = std::min(node.low[d], value);
        node.high[d] = std::max(node.high[d], value);
      }
    }
    int index = static_cast<int>(nodes_.size());
    nodes_.push_back(node);
    if (end - begin <= leaf_size_) return index;
    int split = 0;
    double widest = -1;
    for (int d = 0; d < dimension_; ++d) {
      double span = node.high[d] - node.low[d];
      if (span > widest) { widest = span; split = d; }
    }
    if (!(widest > 0)) return index;
    int middle = begin + (end - begin) / 2;
    std::nth_element(order_.begin() + begin, order_.begin() + middle, order_.begin() + end,
      [&](int a, int b) { double x = at(a, split), y = at(b, split); return x == y ? a < b : x < y; });
    int left = build(begin, middle), right = build(middle, end);
    nodes_[index].left = left; nodes_[index].right = right;
    return index;
  }

  double upper_max_node(int index, const double *q, double best) const {
    const Node &node = nodes_[index];
    if (node.maximum <= best) return best;
    bool fully_inside = true;
    for (int d = 0; d < dimension_; ++d) {
      if (node.high[d] < q[d]) return best;
      if (node.low[d] < q[d]) fully_inside = false;
    }
    if (fully_inside) return std::max(best, node.maximum);
    if (node.left < 0) {
      for (int p = node.begin; p < node.end; ++p) {
        int row = order_[p];
        if (score_[row] <= best) continue;
        bool qualifies = true;
        for (int d = 0; d < dimension_; ++d) if (at(row, d) < q[d]) { qualifies = false; break; }
        if (qualifies) best = score_[row];
      }
      return best;
    }
    int first = node.left, second = node.right;
    if (nodes_[second].maximum > nodes_[first].maximum) std::swap(first, second);
    best = upper_max_node(first, q, best);
    return upper_max_node(second, q, best);
  }

  bool lower_exists_node(int index, const double *q) const {
    const Node &node = nodes_[index];
    bool all_below = true, some_strict = false;
    for (int d = 0; d < dimension_; ++d) {
      if (node.low[d] > q[d]) return false;
      if (node.high[d] > q[d]) all_below = false;
      if (node.low[d] < q[d]) some_strict = true;
    }
    if (all_below && some_strict) return true;
    if (node.left < 0) {
      for (int p = node.begin; p < node.end; ++p) {
        int row = order_[p];
        bool below = true, strict = false;
        for (int d = 0; d < dimension_; ++d) {
          double value = at(row, d);
          if (value > q[d]) { below = false; break; }
          if (value < q[d]) strict = true;
        }
        if (below && strict) return true;
      }
      return false;
    }
    return lower_exists_node(node.left, q) || lower_exists_node(node.right, q);
  }
};

}  // namespace

extern "C" SEXP rd_majorant(SEXP z, SEXP score, SEXP mode, SEXP t, SEXP leaf) {
  const int n = Rf_nrows(z), k = Rf_ncols(z);
  const int leaf_size = Rf_asInteger(leaf);
  Representation rep = represent(REAL(z), n, k, mode_value(mode), Rf_asInteger(t));
  std::map<long long, std::vector<int> > groups;
  for (int i = 0; i < n; ++i) groups[rep.group[i]].push_back(i);
  SEXP answer = PROTECT(Rf_allocVector(REALSXP, n));
  std::copy(REAL(score), REAL(score) + n, REAL(answer));
  for (const auto &entry : groups) {
    if (entry.second.size() < 2) continue;
    DominanceTree tree(rep.coordinate, k, REAL(score), entry.second, leaf_size);
    for (int row : entry.second) {
      REAL(answer)[row] = tree.upper_max(&rep.coordinate[static_cast<size_t>(row) * k], REAL(score)[row]);
    }
  }
  UNPROTECT(1);
  return answer;
}

extern "C" SEXP rd_incon_dec(SEXP z, SEXP rejected, SEXP mode, SEXP t, SEXP leaf) {
  const int n = Rf_nrows(z), k = Rf_ncols(z);
  Representation rep = represent(REAL(z), n, k, mode_value(mode), Rf_asInteger(t));
  std::map<long long, std::vector<int> > selected, omitted;
  for (int i = 0; i < n; ++i) {
    (LOGICAL(rejected)[i] ? selected : omitted)[rep.group[i]].push_back(i);
  }
  SEXP answer = PROTECT(Rf_allocVector(LGLSXP, n));
  std::fill(LOGICAL(answer), LOGICAL(answer) + n, 0);
  for (const auto &entry : selected) {
    auto missing = omitted.find(entry.first);
    if (missing == omitted.end()) continue;
    DominanceTree tree(rep.coordinate, k, nullptr, entry.second, Rf_asInteger(leaf));
    for (int row : missing->second) {
      LOGICAL(answer)[row] = tree.lower_strict_exists(&rep.coordinate[static_cast<size_t>(row) * k]);
    }
  }
  UNPROTECT(1);
  return answer;
}
