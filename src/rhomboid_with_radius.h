#pragma once

// This struct is primarily for testing
template<class FT>
struct RhomboidWithRadius {
  RhomboidWithRadius(Rhomboid rho, FT r) : rho(rho), r(r) {}
  // (Sorted) set of (sorted) vertices.
  Rhomboid rho;
  // Radius/filtration value.
  FT r;
};

// String representation of RhomboidWithRadius
template<class FT>
std::ostream& operator<<(std::ostream& os, RhomboidWithRadius<FT> const& c) {
    os << "(" << c.rho << ", r: " << c.r << ")";
    return os;
}
