#pragma once

// This struct is primarily for testing
template<class FT>
struct CellWithRadius {
  CellWithRadius(std::vector<CVertex> vertices, FT r) : vertices(vertices), r(r) {}
  // (Sorted) set of (sorted) vertices.
  std::vector<CVertex> vertices;
  // Radius/filtration value.
  FT r;
};

// String representation of CellWithRadius
template<class FT>
std::ostream& operator<<(std::ostream& os, CellWithRadius<FT> const& c) {
    os << "(" << c.vertices << ", " << c.r << ")";
    return os;
}
