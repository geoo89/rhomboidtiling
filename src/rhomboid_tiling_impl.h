/*
 * Copyright (c) 2019-2020 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#include <CGAL/Origin.h>

#include "rhomboid_tiling.h"
#include "utils.h"

#include <set>
#include <cassert>
#include <cmath>

/**
 * Constructor for Order-k Delaunay mosaics up to a given order.
 * 
 * Input:
 *    bpoints: vector of input Points.
 *    highest_order: The order up to (including) which to compute the tiling.
 */
// TODO: Refactor into smaller methods.
template<class Dt>
RhomboidTiling<Dt>::RhomboidTiling(
    const std::vector<Point>& bpoints, int highest_order)
        : cur_id(0), highest_order(highest_order), bpoints(bpoints) {
  // squared length of each point when viewed as a vector
  std::vector<FT> squared_lengths;
  for (auto const& p : bpoints) {
    squared_lengths.push_back(Vector(CGAL::ORIGIN, p).squared_length());
  }

  // Step 1: Compute order-1 Delaunay mosaics for k >= 2
  // Make combinatorial vertices of the first-order Delaunay mosaic.
  std::vector<CVertex> vertices = compute_first_order_vertices(bpoints);

  // Turn each point into a weighted point with weight 0.
  std::vector< std::pair<Weighted_point,PIndex> > points;
  for (PIndex i = 0; i < bpoints.size(); ++i) {
    points.push_back(std::make_pair(Weighted_point(bpoints[i], 0), i));
  }

  // Compute first-order Delaunay triangulation as a regular triangulation
  // with weights 0.
  Reg_Tri T(points.begin(), points.end());

  RhomboidMap rhomboids0d = RhomboidMap();
  typename Dt::Regular_triangulation_finite_cells_iterator cit;
  for (cit = Dt::get_finite_cells_begin(T); cit != Dt::get_finite_cells_end(T); ++cit) {
    Rhomboid rho;
    for (int i = 0; i <= dimension; ++i) {
      uint32_t vindex = cit->vertex(i)->info();
      rho.xon.push_back(vertices[vindex][0]);
    }
    std::sort(rho.xon.begin(), rho.xon.end());
    rhomboids0d[rho] = RhomboidInfo<FT>(cur_id++);
  }

  // For dimensions less than d, fill in empty maps for now.
  std::vector<RhomboidMap> rhomboids0 = 
      std::vector<RhomboidMap>(dimension+1, RhomboidMap());
  rhomboids0.push_back(rhomboids0d);
  rhomboids.push_back(rhomboids0);

  // Step 2: Compute order-k Delaunay mosaics for k >= 2
  for (int k = 2; k <= highest_order; ++k) {
    // Set to store the vertex set of the order-k Delaunay mosaic.
    // TODO: Make unordered or sorted vector?
    std::set<CVertex> new_vertices;

    /* Step 2.1: Compute the vertices of the order-k Delaunay mosaic. */
    for (int j = std::max(0, k - dimension); j < k - 1; ++j) {
      for (std::pair<Rhomboid,RhomboidInfo<FT>> element : rhomboids[j][dimension+1]) {
        Rhomboid rho = element.first;
        if (k - j > 1) {
          CCell vxs = rho.get_slice(k - j);
          new_vertices.insert(vxs.begin(), vxs.end());
        }
      }
    }

    // Step 2.2: Compute the first-generation cells of
    // the order-k Delaunay mosaic via a regular triangulation.

    // Store the new vertex set as a vector.
    std::vector<CVertex> new_vertices_vector(new_vertices.begin(), 
                                             new_vertices.end());

    // Step 2.2.1: Construct the geometric vertices and their weights.
    std::vector< std::pair<Weighted_point,PIndex> > new_points;
    int index = 0;
    for (auto const& cv: new_vertices) {
      // We compute the mean of the (lifted) points that are part of the vertex
      // mean: will be the mean of the 3D coordinates.
      Vector mean = CGAL::NULL_VECTOR;
      // mean_sq_length: will be the mean of the 4-th coordinate (i.e. "height")
      FT mean_sq_length = 0;
      for (auto const& pt: cv) {
        mean = mean + Vector(CGAL::ORIGIN, bpoints[pt]);
        mean_sq_length += squared_lengths[pt];
      }
      mean = mean / k;
      mean_sq_length = mean_sq_length / k;
      // mean.squared_length() is the "height" of the mean point in 3D
      // mean_sq_length is the "height" of the mean of the 4D (lifted) points
      // so the weight is how much the mean of the 4D points lies above the
      // paraboloid, because mean.squared_length() is the height of the
      // paraboloid at the mean point in 3D.
      auto weight = mean.squared_length() - mean_sq_length;

      new_points.push_back(std::make_pair(
          Weighted_point(CGAL::ORIGIN + mean, weight),
          index));
      ++index;
    }


    RhomboidMap rhomboidskd = RhomboidMap();

    // Step 2.2.2: Get weighted Delaunay triangulation and indentify its
    // first-generation cells.
    Reg_Tri T(new_points.begin(), new_points.end());
    typename Dt::Regular_triangulation_finite_cells_iterator cit;
    for (cit = Dt::get_finite_cells_begin(T); cit != Dt::get_finite_cells_end(T); ++cit) {
      // TODO: Use pointers here instead of copies?
      CCell vertices;

      for (int i = 0; i < dimension+1; ++i) {
        uint32_t vindex = cit->vertex(i)->info();
        vertices.push_back(new_vertices_vector[vindex]);
      }

      // Check whether the simplex is a first-generation cell.
      CVertex intersec = intersection(vertices);
      if (intersec.size() == k-1) {
        Rhomboid rho;
        rho.xin = intersec;
        CVertex vxunion = setunion(vertices);
        // Compute Xon
        std::set_difference(vxunion.begin(), vxunion.end(),
                            intersec.begin(), intersec.end(),
                            std::inserter(rho.xon, rho.xon.end()));
        std::sort(rho.xon.begin(), rho.xon.end());
        rhomboidskd[rho] = RhomboidInfo<FT>(cur_id++);
      }
    }

    // For dimensions less than d, fill in empty maps for now.
    std::vector<RhomboidMap> rhomboidsk = std::vector<RhomboidMap>(dimension+1, RhomboidMap());
    rhomboidsk.push_back(rhomboidskd);
    rhomboids.push_back(rhomboidsk);
  }

  compute_all_rhomboids();
  compute_all_radii();

}


// Make set of combinatorial vertices of the first-order Delaunay mosaic.
// Each first-order vertex is a singleton set, containing one point index.
template<class Dt>
std::vector<CVertex> RhomboidTiling<Dt>::compute_first_order_vertices(
    const std::vector<Point>& bpoints) {
  std::vector<CVertex> cvertices1;
  for (int i = 0; i < bpoints.size(); ++i) {
    CVertex cv;
    cv.push_back(i);
    cvertices1.push_back(cv);
  }
  return cvertices1;
}


// Compute rhomboids of all dimensions.
// TODO: Think about alternative implementations: 
// - Decrementally go from d to d-1 until 0,
// to omit a lot of duplicates in low dimensions.
// - Implementation to save space by not storing Xin all the time: In
// addition to the above, compute all upper boundaries first, then
// all lower boundaries. Only store reference to Xin of the parent cell
// instead of copy of Xin.
template<class Dt>
void RhomboidTiling<Dt>::compute_all_rhomboids() {
  // Ensure we have maps for lower-dimensional rhomboids of faces at
  // depth > k that we obtain here in the process.
  // TODO: Just don't store (and compute radii) for such faces,
  // except for depth k+1 vertices.
  for (int d = dimension; d >= 0; --d) {
    rhomboids.push_back(std::vector<RhomboidMap>(d+1, RhomboidMap()));
  }

  for (int k = 0; k < highest_order; ++k) {
    for (const std::pair<Rhomboid,RhomboidInfo<FT>>& element : rhomboids[k][dimension+1]) {
      const Rhomboid& rho = element.first;
      assert(rho.dimension() == dimension+1);
      // Each (d+1)-rhomboid has 3^(d+1) faces (including itself).
      for (int fd = 0; fd < std::pow(3, dimension+1); ++fd) {
        Rhomboid face = Rhomboid();
        face.xin = rho.xin;
        face.xon = CVertex();
        int facedim = 0;
        int facedepth = k;
        int fdc = fd;
        // For each point of rho.xon, for the face we can put it in xin, xon or xout.
        for (int i = 0; i < dimension+1; ++i) {
          int ftype = fdc % 3;
          fdc /= 3;
          if (ftype == 0) { // put point in xin
            face.xin.push_back(rho.xon[i]);
            facedepth += 1;
          } else if (ftype == 1) {
            face.xon.push_back(rho.xon[i]);
            facedim += 1;
          } else { // ftype == 2
            //face.xout.push_back(rho.xon[i]); // xout is implicit
          }
        }
        if (facedim != dimension+1) {
          // If facedim == dimension+1, skip because we got face == rho.
          std::sort(face.xin.begin(), face.xin.end());
          std::sort(face.xon.begin(), face.xon.end());
          if (rhomboids[facedepth][facedim].find(face) == rhomboids[facedepth][facedim].end()) {
            // Add it if we don't have it already.
            rhomboids[facedepth][facedim][face] = RhomboidInfo<FT>(cur_id++);
          }
        }
      }
    }
  }
}


template<class Dt>
void RhomboidTiling<Dt>::compute_all_radii() {
  // Set radius for origin to 0
  rhomboids[0][0][Rhomboid()].r = 0;
  // Set radius for edges incident to origin to 0 and depth 1 vertices.
  for (int i = 0; i < (int) bpoints.size(); ++i) {
    Rhomboid rho = Rhomboid();
    rho.xon.push_back(i);
    rhomboids[0][1][rho].r = 0;
    Rhomboid rho2 = Rhomboid();
    rho2.xin.push_back(i);
    rhomboids[1][0][rho2].r = 0;
  }

  for (int d = dimension+1; d >= 2; --d) {
    // Also go to larger k for lower dimensions?
    for (int k = 0; k < highest_order; ++k) {
      for (std::pair<Rhomboid,RhomboidInfo<FT>> element : rhomboids[k][d]) {
        Rhomboid& rho = element.first;
        RhomboidInfo<FT>& info = rhomboids[k][d][rho];
        assert(rho.dimension() == d);
        if (info.r != FT(-1)) continue; // rhomboid is already in an interval

        // Compute radius value of rhomboid.
        Sphere cs = Dt::circumsphere(rho.xon, bpoints);
        info.r = cs.squared_radius();

        // For each of the d Xon points, determine whether it is
        // contained in the vertex that this rhomboid forms an interval with.
        std::vector<bool> is_in_interval_vertex = std::vector<bool>(d, true);
        if (d > 2) {
          // For d == 2 rho is always paired with the bottom vertex.
          // Therefore is_in_interval_vertex is true for all the points
          // and we skip the code below in this case.
          for (int i = 0; i < d; ++i) {
            // Get the circumsphere of Xon with the i-th point removed from it.
            CVertex xon_rest = CVertex(rho.xon);
            xon_rest.erase(xon_rest.begin() + i);
            Sphere cs_new = Dt::circumsphere(xon_rest, bpoints);
            if (cs_new.has_on_bounded_side(bpoints[rho.xon[i]])) {
              // If the i-th point is inside this sphere, then this point is
              // is not part of the interval vertex.
              is_in_interval_vertex[i] = false;
            }
          }
        }

        // rho.xin union the subset of Xon which is_in_interval_vertex
        // forms the vertex that is in the interval with rho.
        // To get all cells in the interval, for each point p from Xon,
        // we can choose whether it stays in Xon, or whether it goes
        // to Xin/Xout (which depends on whether p is_in_interval_vertex).
        for (int fd = 0; fd < std::pow(2, d); ++fd) {
          Rhomboid face = Rhomboid();
          face.xin = rho.xin;
          int fdc = fd;
          for (int i = 0; i < d; ++i) {
            int ftype = fdc % 2;
            fdc /= 2;
            if (ftype == 0) { // put point in xin
              if (is_in_interval_vertex[i]) face.xin.push_back(rho.xon[i]);
              // else it goes to face.xout, which is implicit
            } else { // ftype == 1
              face.xon.push_back(rho.xon[i]);
            }
          }
          if (face.dimension() != d) { // If it's d, skip because we got rho.
            std::sort(face.xin.begin(), face.xin.end());
            std::sort(face.xon.begin(), face.xon.end());
            rhomboids[face.depth()][face.dimension()][face].r = info.r;
          }
        }
      }
    }
  }
}


template<class Dt>
std::vector<Rhomboid>
RhomboidTiling<Dt>::get_rhomboids() {
  std::vector<Rhomboid> rhomboids_out;
  for (int k = 0; k < (int) rhomboids.size(); ++k) {
    if (rhomboids[k].size() > dimension+1) {    
      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& element : rhomboids[k][dimension+1]) {
        const Rhomboid& rho = element.first;
        rhomboids_out.push_back(rho);
      }
    }
  }
  return rhomboids_out;
}


template<class Dt>
std::vector<Rhomboid>
RhomboidTiling<Dt>::get_rhomboids(int order) {
  if (order >= highest_order || order < 0) return std::vector<Rhomboid>();
  if (rhomboids.size() <= order) return std::vector<Rhomboid>();
  if (rhomboids[order].size() <= dimension+1) return std::vector<Rhomboid>();

  std::vector<Rhomboid> rhomboids_out;
  for (const std::pair<Rhomboid,RhomboidInfo<FT>>& element : rhomboids[order][dimension+1]) {
    const Rhomboid& rho = element.first;
    rhomboids_out.push_back(rho);
  }
  return rhomboids_out;
}


template<class Dt>
std::vector<RhomboidWithRadius<typename Dt::FT>>
RhomboidTiling<Dt>::get_rhomboid_filtration() {
  std::vector<RhomboidWithRadius<FT>> rhomboids_out;
  for (int k = 0; k < (int) rhomboids.size(); ++k) {
    for (int d = 0; d < (int) rhomboids[k].size(); ++d) {
      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& element : rhomboids[k][d]) {
        const Rhomboid& rho = element.first;
        const RhomboidInfo<FT>& info = element.second;
        RhomboidWithRadius<FT> rhor(rho, info.r);
        rhomboids_out.push_back(rhor);
      }
    }
  }
  return rhomboids_out;
}


template<class Dt>
std::vector<CCell>
    RhomboidTiling<Dt>::get_slab_mosaic(int order) {
  if (order >= highest_order || order < 0) return std::vector<CCell>();

  std::vector<CCell> cells;
  // (d+1)-cells can be obtained as slabs of (d+1)-rhomboids at d+1 different depths.
  for (int k = order; k >= std::max(0, order-dimension); --k) {
    // For each rhomboid, obtain the slab
    for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][dimension+1]) {
      const Rhomboid& rho = rr.first;
      CCell vertices = rho.get_slab(order - k);
      // Sort by dimension, with lexicographical tiebreaker.
      std::sort(vertices.begin(), vertices.end(),
          [](const CVertex& a, const CVertex& b) -> bool
              { return a.size() != b.size() ? a.size() < b.size() : a < b; });
      cells.push_back(vertices);
    }
  }
  return cells;
}


// The canonical filtration is a list of pairs, each pair consisting
// of a cell and its square radius value. 
// Each cell is a list of combinatorial vertices, with each
// combinatorial vertex being a list of indices into the original
// point set. Each vertex is sorted in ascending order, the vertices in
// each cell are sorted lexicographically, and the whole
// list of cells is sorted lexicographically. Thus the canonical
// representation is unique, and can be used to compare for equality
// when testing the output.
template<class Dt>
std::vector<CellWithRadius<typename Dt::FT>>
    RhomboidTiling<Dt>::get_slice_filtration(int order) {
  if (order > highest_order || order <= 0) return std::vector<CellWithRadius<FT>>();

  std::vector<CellWithRadius<FT>> cells;
  // We treat vertices separately because if we obtain them as slices
  // of edges, we get each vertex multiple times.
  for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[order][0]) {
    const Rhomboid& rho = rr.first;
    const RhomboidInfo<FT>& rho_info = rr.second;
    CCell vertices;
    vertices.push_back(rho.xin);
    cells.push_back(CellWithRadius<FT>(vertices, rho_info.r));
  }
  // Obtain each k-Delaunay d-cell as a slice of a (d+1)-rhomboid
  for (int d = 1; d <= dimension; ++d) {
    // d-cells can be obtained as slices of (d+1)-rhomboids at d+2 different depths.
    // However for k = order and k = order-(d+1) the slice consists of a single vertex,
    // so we omit these.)
    for (int k = order-1; k >= std::max(0, order-d); --k) {
      // For each rhomboid, obtain the slice with its radius value.
      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][d+1]) {
        const Rhomboid& rho = rr.first;
        const RhomboidInfo<FT>& rho_info = rr.second;
        CCell vertices = rho.get_slice(order - k);
        std::sort(vertices.begin(), vertices.end());
        cells.push_back(CellWithRadius<FT>(vertices, rho_info.r));
      }
    }    
  }
  return cells;
}


template<class Dt>
std::vector<CVertex> RhomboidTiling<Dt>::get_vertices(int order) {
  if (order > highest_order || order < 0) return std::vector<CVertex>();

  std::vector<CVertex> vertices;
  for (const std::pair<Rhomboid,RhomboidInfo<FT>>& element : rhomboids[order][0]) {
    const Rhomboid& rho = element.first;
    vertices.push_back(rho.xin);
  }
  return vertices;
}


template<class Dt>
std::vector<CVertex> RhomboidTiling<Dt>::get_vertices() {
  std::vector<CVertex> vertices;
  for (int k = 0; k < highest_order; ++k) {
    std::vector<CVertex> vertices_k = get_vertices(k);
    vertices.insert(vertices.end(), 
        std::make_move_iterator(vertices_k.begin()),
        std::make_move_iterator(vertices_k.end()));
  }
  return vertices;
}


// Returns a list of top-dimensional cells, where
// each cell is a list of combinatorial vertices, with each
// combinatorial vertex being a list of indices into the original
// point set. Each vertex is sorted in ascending order, and the vertices in
// each cell are sorted lexicographically. Thus the output is unique up to
// permutation of the cells.
template<class Dt>
std::vector<CCell>
    RhomboidTiling<Dt>::get_slice_mosaic(int order) {
  if (order > highest_order || order <= 0) return std::vector<CCell>();

  std::vector<CCell> cells;
  for (int k = std::max(0, order - dimension); k < order; ++k) {
    for (std::pair<Rhomboid,RhomboidInfo<FT>> element : rhomboids[k][dimension+1]) {
      Rhomboid rho = element.first;
      CCell vertices = rho.get_slice(order - k);
      std::sort(vertices.begin(), vertices.end());
      cells.push_back(vertices);
    }
  }
  return cells;
}


template<class Dt>
std::vector<CCell>
    RhomboidTiling<Dt>::get_halfint_slice_mosaic(int order) {
  if (order >= highest_order || order < 0) return std::vector<CCell>();

  std::vector<CCell> cells;
  for (int k = std::max(0, order - dimension); k <= order; ++k) {
    for (std::pair<Rhomboid,RhomboidInfo<FT>> element : rhomboids[k][dimension+1]) {
      Rhomboid rho = element.first;
      CCell vertices = rho.get_halfint_slice(order - k);
      std::sort(vertices.begin(), vertices.end());
      cells.push_back(vertices);
    }
  }
  return cells;
}


template<class Dt>
std::vector<CellWithRadius<typename Dt::FT>>
    RhomboidTiling<Dt>::get_slab_filtration(int order) {
  if (order >= highest_order || order < 0) return std::vector<CellWithRadius<FT>>();

  std::vector<CellWithRadius<FT>> cells;
  // Obtain each (d+1)-cell as a slab of a (d+1)-rhomboid
  for (int d = 0; d <= dimension; ++d) {
    // (d+1)-cells can be obtained as slabs of (d+1)-rhomboids at d+1 different depths.
    for (int k = order; k >= std::max(0, order-d); --k) {
      // For each rhomboid, obtain the slab
      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][d+1]) {
        const Rhomboid& rho = rr.first;
        const RhomboidInfo<FT>& rho_info = rr.second;
        CCell vertices = rho.get_slab(order - k);
        // Sort by dimension, with lexicographical tiebreaker.
        std::sort(vertices.begin(), vertices.end(),
            [](const CVertex& a, const CVertex& b) -> bool
                { return a.size() != b.size() ? a.size() < b.size() : a < b; });
        cells.push_back(CellWithRadius<FT>(vertices, rho_info.r));
      }
    }    
  }
  return cells;
}


template<class Dt>
std::vector<BifiltrationCell<typename Dt::FT>>
RhomboidTiling<Dt>::get_bifiltration() {
  return get_bifiltration(1, highest_order);
}


template<class Dt>
std::vector<BifiltrationCell<typename Dt::FT>>
RhomboidTiling<Dt>::get_delaunay_filtration(int order) {
  if (order < 1 || order > highest_order) {
    return std::vector<BifiltrationCell<FT>>();
  }
  return get_bifiltration(order, order);
}


template<class Dt>
std::vector<BifiltrationCell<typename Dt::FT>>
RhomboidTiling<Dt>::get_bifiltration(int minorder, int maxorder) {
  if (minorder < 0) minorder = 0;
  if (maxorder > highest_order) maxorder = highest_order;
  if (minorder > maxorder) return std::vector<BifiltrationCell<FT>>();

  std::vector<BifiltrationCell<FT>> bifiltration;
  for (int k = 0; k <= maxorder; ++k) {
    // Obtain the vertices of the bifiltration
    if (k >= minorder && k <= maxorder) {
      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][0]) {
        const Rhomboid& rho = rr.first;
        const RhomboidInfo<FT>& rho_info = rr.second;
        BifiltrationCell<FT> c;
        c.id = rho_info.id;
        c.d = 0;
        c.r = rho_info.r;
        c.k = k; // depth of the vertex.
        // c.boundary is empty
        bifiltration.push_back(c);
      }
    }
    // The only cells of depth < maxorder that we can obtain from
    // depth-k rhomboids are vertices. Therefore stop here.
    if (k == maxorder) break;

    // Obtain higher dimensional cells of the bifiltration.
    for (int d = 1; d < rhomboids[k].size(); ++d) {
      // Skip if all obtained cells would have depth < minorder
      if (k + d <= minorder) continue;

      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][d]) {
        const Rhomboid& rho = rr.first;
        const RhomboidInfo<FT>& rho_info = rr.second;
        // compute ids of lower boundary cells
        std::vector<int> bdu_ids;
        for (const auto& ru : rho.upper_boundary()) {
          bdu_ids.push_back(rhomboids[k][d-1][ru].id);
        }
        // compute ids of upper boundary cells
        std::vector<int> bdl_ids;
        for (const auto& rl : rho.lower_boundary()) {
          bdl_ids.push_back(rhomboids[k+1][d-1][rl].id);
        }
        // Take slices and slabs of rho. Odd numbers are slabs, even slices.
        // The id of a cell is the id of its parent rhomboid + g*cur_id where
        // g is relative to the rhomboid the cell is a slice/slab of.
        //
        // Illustration below for a 3 dimensional rhomboid with 2 dimensional faces.
        //---------------------------------------------
        //   rho   g  |  upper  g  |  lower  g | depth
        //    o       |    o       |           |  k
        //   / \   1  |   / \   1  |           |  k
        //  o   o  2  |  o   o  2  |      o    |  k+1
        //  |\ /|  3  |   \ /   3  |     /|  1 |  k+1
        //  o o o  4  |    o       |    o o  2 |  k+2
        //   \|/   5  |            |    |/   3 |  k+2
        //    o       |            |    o      |  k+3
        //---------------------------------------------
        for (int g = 1; g < 2*d; ++g) {
          // Skip if the cell we would obtain does not satisfy minorder <= depth <= maxorder.
          if (k + (g+1)/2 > maxorder || k + g/2 < minorder) continue;

          BifiltrationCell<FT> c;
          c.id = g*cur_id + rho_info.id;
          c.r = rho_info.r;
          c.k = k + g/2; // depth of the cell.
          // Slices/slabs of the upper boundary rhomboids.
          if (g < 2*d-2) {
            // The interiors of upper boundary rhomboids are only intersected
            // by the slicing planes up to 2d-3.
            for (int bi : bdu_ids) {
              c.boundary.push_back(g*cur_id + bi);
            }            
          }
          // Slices/slabs of the lower boundary rhomboids.
          if (g > 2) {
            // The interiors of lower boundary rhomboids are only intersected
            // by the slicing planes from 3 onwards.
            for (int bi : bdl_ids) {
              c.boundary.push_back((g-2)*cur_id + bi);
            }
          }
          if (g % 2) {
            // Cell is a slab; their dimension is the same as the sliced rhomboid.
            c.d = d;
            // For slabs, we also need the horizontal boundary cells straight above and below.
            if (d == 1) {
              // For edges, we do want to have the vertices as boundary.
              // However encoding them as a slice of an edge is not a unique
              // representation. Thus we obtain them as vertices (0-dimensional rhomboids),
              // which are the boundary of this edge.
              c.boundary.push_back(bdu_ids[0]);
              c.boundary.push_back(bdl_ids[0]);
            } else {
              // Otherwise, we get these boundary cells as g-1 and g+1 slices of rho.
              // The ifs make sure the slice is a (d-1)-cell and no just a vertex.
              if (g-1 > 0) {
                // boundary cell above
                c.boundary.push_back((g-1)*cur_id + rho_info.id);
              }
              if (g+1 < 2*d) {
                // boundary cell below
                c.boundary.push_back((g+1)*cur_id + rho_info.id);
              }
            }
          } else {
            // Cell is a slice; their dimension is one less than the sliced rhomboid.
            c.d = d-1;
            if (d == 2) {
              // For slices of 2-dimensional rhomboids (i.e. edges), the boundary are
              // slices of 1-dimensional rhomboids (i.e. vertices). These vertices
              // haven't been added yet, because they are not in the interior of those
              // 1-dimensional rhomboids. Instead, we need to find out the actual vertex ids.
              // We get these vertices as the first slice of rho.
              for (const CVertex& vx : rho.get_slice(1)) {
                Rhomboid rho0;
                rho0.xin = vx;
                // Get id of the vertex and add it to the cell's boundary.
                c.boundary.push_back(rhomboids[k+1][0][rho0].id);
              }
            }
          }
          bifiltration.push_back(c);
        } // end for each g
      } // end for each rhomboid
    } // end for each d
  } // end for each k
  return bifiltration;
}



template<class Dt>
std::vector<CCell>
RhomboidTiling<Dt>::get_bifiltration_id_map() {
  // Each rhomboid has up to 2*d slices, with d at most dimension+1.
  // So we need to reserve that much space.
  std::vector<CCell> cells(cur_id*2*(dimension+1));
  for (int k = 0; k < rhomboids.size(); ++k) {
    // Obtain the vertices of the bifiltration
    for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][0]) {
      const Rhomboid& rho = rr.first;
      const RhomboidInfo<FT>& rho_info = rr.second;
      cells[rho_info.id] = rho.get_slice(0);
    }
    // Obtain higher dimensional cells of the bifiltration.
    for (int d = 1; d < rhomboids[k].size(); ++d) {
      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][d]) {
        const Rhomboid& rho = rr.first;
        const RhomboidInfo<FT>& rho_info = rr.second;
        for (int g = 1; g < 2*d; ++g) {
          if (g % 2) {
            cells[rho_info.id + cur_id*g] = rho.get_slab(g / 2);
          } else {
            cells[rho_info.id + cur_id*g] = rho.get_slice(g / 2);
          }
        }
      }
    }
  }
  // Make each cell canonical.
  for (auto& cell : cells) {
    std::sort(cell.begin(), cell.end());
  }
  return cells;
}


template<class Dt>
std::vector<BifiltrationCell<typename Dt::FT>>
RhomboidTiling<Dt>::get_unsliced_bifiltration() {
  return get_unsliced_bifiltration(highest_order);
}


// TODO: refactor common functionality with get_bifiltration?
template<class Dt>
std::vector<BifiltrationCell<typename Dt::FT>>
RhomboidTiling<Dt>::get_unsliced_bifiltration(int maxorder) {
  if (maxorder > highest_order) maxorder = highest_order;
  if (maxorder < 0) return std::vector<BifiltrationCell<FT>>();

  std::vector<BifiltrationCell<FT>> bifiltration;
  // Iterate over depth (of the anchor vertex of the rhomboid) k.
  for (int k = 0; k <= maxorder; ++k) {
    // Iterate over dimension d.
    for (int d = 0; d < rhomboids[k].size(); ++d) {
      // For rhomboids anchored at depth maxorder, their slice at maxorder
      // can only be a vertex. But those vertices are already represented by
      // 0-dimensional rhomboids.
      if ((k == maxorder) && (d != 0)) break;
      // Collect all rhomboids of depth k and dimension d.
      for (const std::pair<Rhomboid,RhomboidInfo<FT>>& rr : rhomboids[k][d]) {
        const Rhomboid& rho = rr.first;
        const RhomboidInfo<FT>& rho_info = rr.second;
        BifiltrationCell<FT> c;
        c.id = rho_info.id;
        c.r = rho_info.r;
        c.d = d;
        c.k = k;

        // compute ids of upper boundary cells
        std::vector<int> bdu_ids;
        for (const auto& ru : rho.upper_boundary()) {
          bdu_ids.push_back(rhomboids[k][d-1][ru].id);
          c.boundary.push_back(rhomboids[k][d-1][ru].id);
        }

        std::vector<int> bdl_ids;
        if ((k + 1 < maxorder) || (d == 1)) {
          // Lower boundary cells are of depth k+1.
          // So if k + 1 == maxorder, we discard them,
          // unless their are vertices (d == 1).
          // Compute ids of lower boundary cells.
          for (const auto& rl : rho.lower_boundary()) {
            bdl_ids.push_back(rhomboids[k+1][d-1][rl].id);
            c.boundary.push_back(rhomboids[k+1][d-1][rl].id);
          }
        }

        if (k + d > maxorder) {
          // This rhomboid got sliced at maxorder.
          // Note that this can only happen for rhomboids of d >= 2.
          // We add its lower boundary, encoded with the following ID:
          c.boundary.push_back(rho_info.id + cur_id);
          // This is a new cell, so we also need to compute its stats...
          BifiltrationCell<FT> c2;
          c2.id = rho_info.id + cur_id;
          c2.r = rho_info.r;
          c2.d = d-1;
          c2.k = maxorder;
          // ...and its boundary.
          int g = 2*(maxorder - k);
          // Illustration below for a 3 dimensional rhomboid with 2 dimensional faces.
          //---------------------------------------------
          //   rho   g  |  upper  g  |  lower  g | depth
          //    o       |    o       |           |  k
          //   / \   1  |   / \   1  |           |  k
          //  o   o  2  |  o   o  2  |      o    |  k+1
          //  |\ /|  3  |   \ /   3  |     /|  1 |  k+1
          //  o o o  4  |    o       |    o o  2 |  k+2
          //   \|/   5  |            |    |/   3 |  k+2
          //    o       |            |    o      |  k+3
          //---------------------------------------------
          if (g < 2*d-2) {
            // Slices of the upper boundary rhomboids.
            // The interiors of upper boundary rhomboids are only intersected
            // by the slicing planes up to 2d-3.
            for (int bi : bdu_ids) {
              c2.boundary.push_back(cur_id + bi);
            }            
          }
          if (g > 2) {
            // Slices of the lower boundary rhomboids.
            // The interiors of lower boundary rhomboids are only intersected
            // by the slicing planes from 3 onwards.
            for (int bi : bdl_ids) {
              c2.boundary.push_back(cur_id + bi);
            }
          }
          if (d == 2) {
            // For slices of 2-dimensional rhomboids (i.e. edges), the boundary are
            // slices of 1-dimensional rhomboids (i.e. vertices). These vertices
            // haven't been added yet, because they are not in the interior of those
            // 1-dimensional rhomboids. Instead, we need to find out the actual vertex ids.
            // We get these vertices as the first slice of rho.
            for (const CVertex& vx : rho.get_slice(1)) {
              Rhomboid rho0;
              rho0.xin = vx;
              // Get id of the vertex and add it to the cell's boundary.
              c2.boundary.push_back(rhomboids[maxorder][0][rho0].id);
            }
          }
          bifiltration.push_back(c2);
        }

        bifiltration.push_back(c);
      } // end for each rhomboid
    } // end for each d
  } // end for each k
  return bifiltration;
}

