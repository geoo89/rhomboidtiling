#ifndef _TYPES_H_
#define _TYPES_H_

#include <vector>

// point index
typedef unsigned            PIndex;
// A combinatorial vertex is a set of point indices,
// represented as a sorted vector of point indices.
typedef std::vector<PIndex>  CVertex;
// A combinatorial cell is a set of CVertices,
// represented as a sorted vector of CVertices.
typedef std::vector<CVertex>  CCell;

#endif // _TYPES_H_
