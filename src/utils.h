#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "types.h"

// Compute common intersection of a set of combinatorial vertices, which
// are sets of point indices.
CVertex intersection(const std::vector<CVertex> &vertices);

// Compute common intersection of a set of combinatorial vertices, which
// are sets of point indices.
CVertex setunion(const std::vector<CVertex> &vertices);

// Subsets of size j for set on i elements. More specifically,
// combinatorial_subsets[i][j] contains a vector of all subsets of {0, ..., i-1}
// of size j, for i up to 4 and j <= i. Each subset is represented as a sorted vector.
// Used to get second and third generation cells from a first generation cell,
// and for getting Delaunay slices from rhomboids.
extern const std::vector<std::vector<std::vector<std::vector<int> > > >
    combinatorial_subsets;

// Give nice string representation of vectors when streaming to output. 
template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  if (!v.empty()) {
    out << '[';
    for (unsigned i = 0; i < v.size() - 1; ++i) {
      out << v[i] << ", ";
    }
    out << v[v.size() - 1] << "]";
  } else {
    out << "[]";
  }

  return out;
}

namespace std {
  // Hash function for vectors.
  template <typename T>
  struct hash<std::vector<T>> {
    std::size_t operator()(const std::vector<T>& v) const {
      std::size_t seed = v.size();
      for (const T& i : v) 
        seed ^= std::hash<T>()(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
      return seed;
    }
  };

  // Based on the boost implementation. seed is modified inplace.
  template<typename T>
  void hash_combine(size_t& seed, const T& value) {
    seed ^= std::hash<T>()(value) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  };

} // namespace std
