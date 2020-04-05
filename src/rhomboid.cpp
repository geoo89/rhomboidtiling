#include "rhomboid.h"
#include "utils.h"

std::vector<CVertex> Rhomboid::get_slice(int depth) const {
  std::vector<CVertex> vertices;
  int d = xon.size();
  // Xin, together with any depth out of the d+1 Xon points, make a vertex.
  for (const auto& subset : combinatorial_subsets[d][depth]) {
    CVertex vx = CVertex(xin);
    for (const auto& element : subset) {
      // add element to vertex
      vx.push_back(xon[element]);
    }
    std::sort(vx.begin(), vx.end());
    vertices.push_back(vx);
  }
  return vertices;
}

CCell Rhomboid::get_slab(int depth) const {
  CCell vertices = get_slice(depth);
  CCell vertices2 = get_slice(depth+1);
  // concatenate the depth and depth+1 vertex sets.
  vertices.insert(vertices.end(), 
      std::make_move_iterator(vertices2.begin()),
      std::make_move_iterator(vertices2.end()));
  return vertices;
}

CCell Rhomboid::get_halfint_slice(int depth) const {
  CCell vertices;
  CCell edge_bottoms = get_slice(depth+1);
  for (const auto& vertex : edge_bottoms) {
    for (int i = 0; i < (int) vertex.size(); ++i) {
      CVertex vertex_copy = CVertex(vertex);
      // remove an element to get the top vertex of the edge
      vertex_copy.erase(vertex_copy.begin() + i);
      // add the removed point index back at the end to make the bottom vertex
      vertex_copy.push_back(vertex[i]);
      vertices.push_back(vertex_copy);
    }
  }
  return vertices;
}

// Those boundary rhomboids that share the same anchor vertex.
// Same Xin, and Xon minus any vertex.
std::vector<Rhomboid> Rhomboid::upper_boundary() const {
  std::vector<Rhomboid> bd;
  for (int i = 0; i < dimension(); ++i) {
    Rhomboid rho;
    rho.xin = xin;
    rho.xon = xon;
    rho.xon.erase(rho.xon.begin()+i);
    bd.push_back(rho);
  }
  return bd;
}

// Those boundary rhomboids that do not share the same anchor vertex.
// Their anchor vertex is one level deeper.
std::vector<Rhomboid> Rhomboid::lower_boundary() const {
  std::vector<Rhomboid> bd;
  for (int i = 0; i < dimension(); ++i) {
    Rhomboid rho;
    rho.xin = xin;
    rho.xon = xon;
    PIndex pi = xon[i];
    rho.xon.erase(rho.xon.begin()+i);
    rho.xin.push_back(pi);
    std::sort(rho.xin.begin(), rho.xin.end());
    bd.push_back(rho);
  }
  return bd;
}

// Xon minus any vertex, Xin union that vertex.
std::vector<Rhomboid> Rhomboid::boundary() const {
  std::vector<Rhomboid> bd = upper_boundary();
  std::vector<Rhomboid> bd2 = lower_boundary();
  // concatenate the depth and depth+1 vertex sets.
  bd.insert(bd.end(), 
      std::make_move_iterator(bd2.begin()),
      std::make_move_iterator(bd2.end()));
  return bd;
}

// String representation of Rhomboid
std::ostream& operator<<(std::ostream& os, const Rhomboid& rho) {
    os << "Xin: " << rho.xin << ", Xon: " << rho.xon;
    return os;
}
