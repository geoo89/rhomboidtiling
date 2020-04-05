#ifndef _RHOMBOID_H_
#define _RHOMBOID_H_

#include "types.h"
#include "utils.h"

// A Rhomboid.
// The dimension is implicit as xon.size().
// The depth is implicit as xin.size().
struct Rhomboid {
  // Indices of the points that all vertices in this rhomboid share.
  CVertex xin;
  // Indices of the points that are present in some vertices of this rhomboid.
  CVertex xon;
  // Sufficient subset of Xout fully defining the sphere of the rhomboid.
  //CVertex xout;

  Rhomboid() : xin(CVertex()), xon(CVertex()) {}
  Rhomboid(CVertex xin, CVertex xon) : xin(xin), xon(xon) {}

  int dimension() const { return (int) xon.size(); }
  int depth() const { return (int) xin.size(); }

  // Get the slice at relative depth, which is a higher order Delaunay cell.
  CCell get_slice(int depth) const;
  // Get the slab at relative depth.
  CCell get_slab(int depth) const;
  // Get the half-integers slice at relative depth + 0.5,
  // which is dual to a degree-k Voronoi cell.
  // Its vertices are represented as CVertex (vectors of point indices),
  // with a vertex being the slice of the edge between all but the last
  // element of the CVertex and the full CVertex.
  // Example: [0, 2, 3, 1] is the edge going from [0, 2, 3] to [0, 1, 2, 3].
  CCell get_halfint_slice(int depth) const;

  // Vector of boundary rhomboids.
  std::vector<Rhomboid> boundary() const;
  // Those boundary rhomboids that share the same anchor vertex.
  std::vector<Rhomboid> upper_boundary() const;
  // Those boundary rhomboids that do not share the same anchor vertex.
  // Their anchor vertex is one level deeper.
  std::vector<Rhomboid> lower_boundary() const;
};

// Needed if using a map rather than unordered_map of rhomboids.
inline int operator<(const Rhomboid& r0, const Rhomboid& r1) {
    return std::make_pair(r0.xin, r0.xon) < std::make_pair(r1.xin, r1.xon);
}

// Needed for unordered_map of rhomboids.
inline int operator==(const Rhomboid& r0, const Rhomboid& r1) {
    return r0.xin == r1.xin && r0.xon == r1.xon;
}

// String representation of Rhomboid
std::ostream& operator<<(std::ostream& os, const Rhomboid& rho);

// Hash function for rhomboids.
// Needed for unordered_map of rhomboids.
namespace std {
  template <>
  struct hash<Rhomboid> {
    std::size_t operator()(const Rhomboid& rho) const {
      std::size_t seed = std::hash<CVertex>()(rho.xin);
      hash_combine(seed, rho.xon);
      return seed;
    }
  };
} // namespace std

template<class FT>
struct RhomboidInfo {
  RhomboidInfo() : r(FT(-1)), id(-1) {}
  RhomboidInfo(int id) : r(FT(-1)), id(id) {}
  // RhomboidInfo(CVertex xo, double r, int id) : xo(xo), r(r), id(id) {}
  // The interval vertex is the vertex that has the same radius value as
  // the rhomboid. It is obtained as xin union xo.
  //CVertex xo;
  // Radius/filtration value.
  FT r;
  // An id unique across all rhomboids.
  int id;
};

#endif // _RHOMBOID_H_
