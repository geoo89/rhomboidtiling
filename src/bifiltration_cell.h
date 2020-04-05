#pragma once

// Class representing a cell in the Delaunay bifiltration.
// Each cell has an id, filtration values and ids of its boundary cells.
template<class FT>
struct BifiltrationCell {
  // Unique ID of the cell.
  // This is the ID of the rhomboid the cell is a slice of,
  // plus 2 times the relative depth of the slice times the total number of
  // rhomboids (for this purpose, we consider the depth of a slab to be the
  // mean of the depths of its two bounding slices):
  //   o   0
  //  / \  0.5
  // o   o 1
  //  \ /  1.5
  //   o   2
  // Vertices are obtained as vertices rather than slices of edges, to ensure
  // uniqueness.
  int id;
  // Radius filtration value.
  FT r;
  // Depth of the cell.
  int k;
  // Dimension of the cell.
  int d;
  // ids of the boundary cells.
  std::vector<int> boundary;
};

// String representation of BifiltrationCell
template<class FT>
std::ostream& operator<<(std::ostream& os, BifiltrationCell<FT> const& c) {
    os << c.id << " " << c.d << " " << c.k << " " << c.r;
    for (int bd : c.boundary) {
      os << " " << bd;
    }
    return os;
}

// For sorting convenience. Dimension first, then k, then r, then id.
template<class FT>
int operator<(const BifiltrationCell<FT>& r0, const BifiltrationCell<FT>& r1) {
    return r0.d != r1.d ? r0.d < r1.d :
           r0.k != r1.k ? r0.k < r1.k :
           r0.r != r1.r ? r0.r < r1.r :
           r0.id < r1.id;
}
