#pragma once

// Like BifiltrationCell, but cell ids are replaced combinatorial cells.
// This is to ensure uniqueness to allow comparison.
template<class FT>
struct CombinatorialBifiltrationCell {
    CombinatorialBifiltrationCell(int d, int k, std::vector<CVertex> vertices,
                            FT r, std::vector<std::vector<CVertex>> boundary) :
            d(d), k(k), vertices(vertices), r(r), boundary(boundary) { }

    CombinatorialBifiltrationCell(BifiltrationCell<FT> bc, std::vector<std::vector<CVertex>> cell_id_map) :
            d(bc.d), k(bc.k), r(bc.r) {
        vertices = cell_id_map[bc.id];
        for (const auto& bd : bc.boundary) {
            boundary.push_back(cell_id_map[bd]);
        }
    }

    // combinatorial cell (set of vertices)
    std::vector<CVertex> vertices;
    // Dimension of the cell.
    int d;
    // Depth of the cell.
    int k;
    // Filtration value.
    FT r;
    // combinatorial boundary cells.
    std::vector<std::vector<CVertex>> boundary;
};


// To define unique ordering.
template<class FT>
int operator<(const CombinatorialBifiltrationCell<FT>& r0, const CombinatorialBifiltrationCell<FT>& r1) {
    return r0.d != r1.d ? r0.d < r1.d :
           r0.k != r1.k ? r0.k < r1.k :
           r0.vertices < r1.vertices;
}

// String representation of CombinatorialBifiltrationCell
template<class FT>
std::ostream& operator<<(std::ostream& out, const CombinatorialBifiltrationCell<FT>& c) {
  out << "d:" << c.d 
      << ", k:" << c.k 
      << ", vxs: " << c.vertices 
      << ", r:" << c.r 
      << ", bd:" << c.boundary;
  return out;
}
