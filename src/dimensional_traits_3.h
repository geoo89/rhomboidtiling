#pragma once

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
// CGAL <4.10
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <vector>

#include "types.h"


template<class K>
class DimensionalTraits_3 {
  public:
    static const int dimension = 3;

    typedef typename K::Point_3                                                       Point;
    typedef typename K::Vector_3                                                     Vector;
    typedef typename K::Sphere_3                                                     Sphere;
    typedef typename K::FT                                                               FT;

    // CGAL 4.10 breaks some backward compatibility.
    // This code is for CGAL <4.10, but replacing the typedefs below
    // should make to code compile for later CGAL versions.

    // If using CGAL <4.10:
    typedef typename CGAL::Regular_triangulation_euclidean_traits_3<K>               Traits;
    typedef typename CGAL::Triangulation_vertex_base_3<Traits>                        Vbase;
    typedef typename CGAL::Triangulation_vertex_base_with_info_3<PIndex, Traits,Vbase>   Vb;
    typedef typename CGAL::Regular_triangulation_cell_base_3<Traits>                     Cb;
    typedef typename CGAL::Triangulation_data_structure_3<Vb,Cb>                        Tds;
    typedef typename CGAL::Regular_triangulation_3<Traits, Tds>       Regular_triangulation;
    typedef typename Traits::Weighted_point_3                                Weighted_point;
    typedef typename Regular_triangulation::Finite_cells_iterator
                                                Regular_triangulation_finite_cells_iterator;

    // If using CGAL >=4.10:
    // typedef typename K::FT                                               Weight;
    // typedef typename K::Weighted_point_3                                 Weighted_point;
    // typedef typename CGAL::Regular_triangulation_vertex_base_3<K>        Vb0;
    // typedef typename CGAL::Triangulation_vertex_base_with_info_3<PIndex, K, Vb0> Vb;
    // typedef typename CGAL::Regular_triangulation_cell_base_3<K>          Cb;
    // typedef typename CGAL::Triangulation_data_structure_3<Vb,Cb>         Tds;
    // typedef typename CGAL::Regular_triangulation_3<K, Tds>               Reg_Tri;

    static Regular_triangulation_finite_cells_iterator
    get_finite_cells_begin(const Regular_triangulation& T) {return T.finite_cells_begin();}
    static Regular_triangulation_finite_cells_iterator
    get_finite_cells_end(const Regular_triangulation& T) {return T.finite_cells_end();}

    // Get circumsphere of the points with the given indices.
    static Sphere
    circumsphere(const CVertex& pids, const std::vector<Point>& bpoints) {
        if (pids.size() == 4) {
            return Sphere(bpoints[pids[0]], bpoints[pids[1]], bpoints[pids[2]], bpoints[pids[3]]);
        } else if (pids.size() == 3) {
            return Sphere(bpoints[pids[0]], bpoints[pids[1]], bpoints[pids[2]]); 
        } else if (pids.size() == 2) {
            return Sphere(bpoints[pids[0]], bpoints[pids[1]]); 
        } else {
            // return std::optional<Sphere>() or an exception in the future.
            // Currently there are no cases where this should happen.
            return Sphere(Point(CGAL::ORIGIN));
        }
    }

    // Make a point with the given coordinates.
    static Point
    make_point(const std::vector<double>& coords) {
        if (coords.size() == 3) {
            return Point(coords[0], coords[1], coords[2]);
        } else {
            // return std::optional<Point>() or an exception in the future.
            return Point(CGAL::ORIGIN);
        }
    }
};
