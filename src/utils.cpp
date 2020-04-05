
#include "utils.h"

// Compute common intersection of a set of combinatorial vertices, which
// are sets of point indices.
CVertex intersection(const std::vector<CVertex> &vertices) {
    CVertex previous_intersection = vertices[0];
    CVertex new_intersection;

    for (unsigned i = 1; i < vertices.size(); ++i) {
        std::set_intersection(
            previous_intersection.begin(), previous_intersection.end(),
            vertices[i].begin(), vertices[i].end(),
            std::back_inserter(new_intersection));
        std::swap(previous_intersection, new_intersection);
        new_intersection.clear();
    }
    return previous_intersection;
}


// Compute common intersection of a set of combinatorial vertices, which
// are sets of point indices.
CVertex setunion(const std::vector<CVertex> &vertices) {
    CVertex previous_union = vertices[0];
    CVertex new_union;

    for (unsigned i = 1; i < vertices.size(); ++i) {
        std::set_union(
            previous_union.begin(), previous_union.end(),
            vertices[i].begin(), vertices[i].end(),
            std::back_inserter(new_union));
        std::swap(previous_union, new_union);
        new_union.clear();
    }
    return previous_union;
}


// All subsets of {0, ..., i} for i up to 3, indexed by their size.
// Used to get second and third generation cells from a first generation cell,
// and for getting Delaunay slices from rhomboids.
// Ordering of the subsets is lexicographic.
const std::vector<std::vector<std::vector<int> > >
    combinatorial_subsets4 = {
      {{}}, // the trivial subset
      {{0},{1},{2},{3}}, 
      {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}, 
      {{0,1,2},{0,1,3},{0,2,3},{1,2,3}},
      {{0,1,2,3}}
    };
const std::vector<std::vector<std::vector<int> > >
    combinatorial_subsets3 = {
      {{}},
      {{0},{1},{2}}, 
      {{0,1},{0,2},{1,2}},
      {{0,1,2}}
    };
const std::vector<std::vector<std::vector<int> > >
    combinatorial_subsets2 = {
      {{}},
      {{0},{1}}, 
      {{0,1}}
    };
const std::vector<std::vector<std::vector<int> > >
    combinatorial_subsets1 = {
      {{}},
      {{0}}
    };

const std::vector<std::vector<std::vector<std::vector<int> > > >
    combinatorial_subsets = {
      {{{}}},
      combinatorial_subsets1,
      combinatorial_subsets2,
      combinatorial_subsets3,
      combinatorial_subsets4
    };
