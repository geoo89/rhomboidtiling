#pragma once

#include "catch2/catch.hpp"

using namespace Catch;

// The following are defined in here instead of the respective header files
// because we use Catch2's Approx functionality to compare the radius.
// TODO: We shouldn't need Approx on the radius.
// Figure out why the values are not exact.

// Comparison operator for CellWithRadius, for Catch2
template<class FT>
bool operator==(const CellWithRadius<FT>& c1, const CellWithRadius<FT>& c2) {
    return c1.vertices == c2.vertices && CGAL::to_double(c1.r) == Approx(CGAL::to_double(c2.r));
}

// Comparison operator for RhomboidWithRadius, for Catch2
template<class FT>
bool operator==(const RhomboidWithRadius<FT>& c1, const RhomboidWithRadius<FT>& c2) {
    return c1.rho == c2.rho && CGAL::to_double(c1.r) == Approx(CGAL::to_double(c2.r));
}

// Comparison operator for CombinatorialBifiltrationCell, for Catch2
template<class FT>
int operator==(const CombinatorialBifiltrationCell<FT>& r0, const CombinatorialBifiltrationCell<FT>& r1) {
    return r0.d == r1.d && r0.k == r1.k && r0.vertices == r1.vertices
           && r0.boundary == r1.boundary && CGAL::to_double(r0.r) == Approx(CGAL::to_double(r1.r));
}
