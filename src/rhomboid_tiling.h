/*
 * Copyright (c) 2019-2020 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#ifndef _RHOMBOID_TILING_H_
#define _RHOMBOID_TILING_H_

#include <vector>
#include <unordered_map>

#include "rhomboid.h"
#include "cell_with_radius.h"
#include "rhomboid_with_radius.h"
#include "bifiltration_cell.h"

template<class Dt>
class RhomboidTiling {
  public:
    typedef typename Dt::Point                                                   Point;
    typedef typename Dt::Vector                                                 Vector;
    typedef typename Dt::Sphere                                                 Sphere;
    typedef typename Dt::FT                                                         FT;
    typedef typename Dt::Regular_triangulation                                 Reg_Tri;
    typedef typename Dt::Weighted_point                                 Weighted_point;

    // To save memory, use sorted vector instead, with binary search to find boundaries?
    typedef typename std::unordered_map<Rhomboid,RhomboidInfo<FT>>         RhomboidMap;

    static const int dimension = Dt::dimension;

    /**
     * Constructor for RhomboidTiling up to a given order (depth).
     *
     * Upon creation, this class computes the tiling up to the specified depth,
     * including rhomboids of all dimensions and their radius values.
     * 
     * Input:
     *    bpoints: vector of input Points.
     *    highest_order: The order up to which to compute the
     *        rhomboid tiling. This means all rhomboids that have a non-trivial
     *        slice at highest_order are included, so in particular all
     *        rhomboids whose anchor vertex is at depth highest_order-1.
     *        Thus in particular the
     *        order-k Delaunay mosaics up to including highest_order can be
     *        obtained from the tiling.
     */
    RhomboidTiling(const std::vector<Point>& bpoints, int highest_order);
    /* Return the top-dimensional cells of the rhomboid tiling.
     *
     * Each cell is a Rhomboid object, consisting of a set of point indices
     * Xin and a set Xon, with the set Xout being implicit.
     */
    std::vector<Rhomboid> get_rhomboids();
    /* Return the top-dimensional cells of the rhomboid tiling whose anchor
     * vertex is at the specified depth (order).
     */
    std::vector<Rhomboid> get_rhomboids(int order);
    /* For the specified order k, return the set of combinatorial vertices.
     *
     * Each vertex is represented as a k-tuple of point indices.
     *
     * Its geometric location can be obtained as the barycenter of the
     * k points that the combinatorial vertex refers to.
     */
    std::vector<CVertex> get_vertices(int order);
    // Get all vertices up to (including) highest_order.
    std::vector<CVertex> get_vertices();
    /* For the specified order k, get the top-dimensional cells of the
     * order-k Delaunay mosaic.
     *
     * Each cell is represented in as a sorted tuple of combinatorial vertices,
     * with each combinatorial vertex being a sorted tuple of k indices into
     * the original vector of input points. The returned vector is unique up to
     * permutation.
     *
     * Thus two order-k Delaunay mosaics are combinatorially equivalent if and
     * only if their sorted canonical representations are equivalent.
     * This is useful for testing.
     */
    std::vector<CCell> get_slice_mosaic(int order);
    /* For the specified order k, get the top-dimensional cells of the slice of
     * the rhomboid tiling at depth k + 0.5, which is dual to the degree-k
     * Voronoi tesselation.
     *
     * Each cell is represented in as a sorted tuple of combinatorial vertices,
     * with each combinatorial vertex being a vector of k+1 point indices,
     * encoding the endpoint of the edge the vertex is a slice of.
     * The the first k elements of the vector are in sorted order and encode the
     * top vertex of the edge. The last point index is the point index to be
     * added to the top vertex to obtain the bottom vertex.
     * The returned vector is unique up to permutation.
     */
    std::vector<CCell> get_halfint_slice_mosaic(int order);
    /* For the specified order k, get the canonical representation (up to
     * permutation of cells) of top-dimensional cells of the slab between
     * the order-k Delaunay mosaic and the order-(k+1) Delaunay mosaic.
     */
    std::vector<CCell> get_slab_mosaic(int order);

    /* Return the filtration of the rhomboid tiling.
     * 
     * This returns all rhomboids (of any dimensions) with their respective
     * filtration value (radius r).
     */
    std::vector<RhomboidWithRadius<FT>> get_rhomboid_filtration();
    /* For the specified order k, get the canonical representation of the
     * order-k Delaunay mosaic with radius values (cells of all dimensions).
     */
    std::vector<CellWithRadius<FT>> get_slice_filtration(int order);
    /* For the specified order k, get the canonical representation of the slab
     * between the order-k Delaunay mosaic and the order-(k+1) Delaunay mosaic
     * with radius values (cells of all dimensions). Only cells exclusive to
     * the slab are returned. In particular that means the result does not
     * contain any cells from the order-k Delaunay mosaic and the order-(k+1)
     * Delaunay mosaic, it contains no vertices and it is not closed under
     * taking faces.
     */
    std::vector<CellWithRadius<FT>> get_slab_filtration(int order);

    // Get the rhomboid bifiltration with parameters r (radius) and
    // k (depth of the anchor vertex of the rhomboid) up to highest_order
    // (for cells of all dimensions). Each cell is a BifiltrationCell<FT>.
    // It is assigned an ID and has radius, depth,
    // dimension and the ids of its boundary cells associated to it.
    std::vector<BifiltrationCell<FT>> get_rhomboid_bifiltration();
    // The rhomboid bifiltration between minorder and maxorder (inclusive).
    std::vector<BifiltrationCell<FT>> get_rhomboid_bifiltration(int minorder, int maxorder);
    // Get the bifiltration with parameters r (radius) and k (depth)
    // up to highest_order (cells of all dimensions). Each cell is a
    // BifiltrationCell<FT>. It is assigned an ID and has radius, depth,
    // dimension and the ids of its boundary cells associated to it.
    // To get the combinatorial representation of a cell, use the method
    // get_bifiltration_id_map().
    std::vector<BifiltrationCell<FT>> get_bifiltration();
    // The bifiltration between minorder and maxorder (inclusive).
    std::vector<BifiltrationCell<FT>> get_bifiltration(int minorder, int maxorder);
    // Get order-k Delaunay filtration for the specified order k.
    // This is the same as get_slice_filtration, however with a different
    // output format convenient for creating a boundary matrix.
    std::vector<BifiltrationCell<FT>> get_delaunay_filtration(int order);
    // A vector of cells of the bifiltration.
    // The entry at position j is the cell with id j.
    // As the cell ids are not continuous, many entries of the vector
    // will be blank and not refer to a cell.
    std::vector<CCell> get_bifiltration_id_map();

  private:
    // Make set of combinatorial vertices of the first-order Delaunay mosaic.
    // Each first-order vertex is a singleton set, containing one point index.
    std::vector<CVertex> compute_first_order_vertices(const std::vector<Point>& bpoints);
    // TODO: Evaluate these lazily?
    // From the top dimensional rhomboids, compute all lower dimensional ones.
    void compute_all_rhomboids();
    // Compute filtration values for all rhomboids.
    void compute_all_radii();

    // The set of points to compute the rhomboid tiling of
    std::vector<Point> bpoints;
    // id counter to assign unique ids to all rhomboids.
    int cur_id;
    // Up to which order to compute order-k Delaunay. Equivalently, all
    // rhomboids with anchor up to depths highest_order-1 are computed.
    int highest_order;
    // A 2D array (indexed by k and d) of unordered_maps mapping each rhomboid
    // to additional data.
    std::vector<std::vector<RhomboidMap> > rhomboids;
};

#include "rhomboid_tiling_impl.h"

#endif // _RHOMBOID_TILING_H_
