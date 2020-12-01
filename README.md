# Rhomboid Tiling and Order-k Delaunay Mosaics

For a given d-dimensional finite set of input points, the rhomboid tiling,
introduced in [The Multi-cover Persistence of Euclidean
Balls](http://pub.ist.ac.at/~edels/Papers/2018-P-01-MultiCover.pdf),
is a (d+1)-dimensional geometric structure consisting of rhomboids, which
are parallelpipeds and their analogues in dimensions other than 3.
From the rhomboid tiling one can obtain order-k Delaunay mosaics (dual to
order-k Voronoi tesselations), degree-k Delaunay mosaics (dual to degree-k
Voronoi tesselations, which decompose the space into k-th Brillouin zones),
and a discrete version of the filtration of multi-covers of Euclidean balls.

This is a C++ implementation of an algorithm to compute the rhomboid tiling
for a given input point set, based on the paper 
[A Simple Algorithm for Computing Higher Order Delaunay
Mosaics](http://pub.ist.ac.at/~edels/Papers/2020-P-01-SimpleAlgorithm.pdf).
Currently, the input point set can be 2- or 3-dimensional.
The implementation uses the [CGAL](https://www.cgal.org/)
library for Delaunay triangulations and exact arithmetics and includes
some unit tests using [Catch2](https://github.com/catchorg/Catch2),
and mostly makes an earlier implementation from 
https://github.com/geoo89/orderkdelaunay obsolete.


## Prerequisites

_Prerequisites:_ cmake, CGAL version <= 4.9, Catch2 (included);
to work with CGAL version >= 4.10, some typedefs need to be changed,
see `src/dimensional_traits_2.h` and `src/dimensional_traits_3.h`.

The build setup builds a commandline tool and tests. To build, run:
```
cmake .
make
```

The unit tests use Catch2 which is included as a header and comes
with its own licence.


## Interface

The class representing the rhomboid tiling is RhomboidTiling<Dt> as defined in
`src/rhomboid_tiling.h`. It is templated with a dimensional traits class that
contains dimension specific functionality. For 2D computations use 
DimensionalTraits_2 from `src/dimensional_traits_2.h`; for 3D use
DimensionalTraits_3 from `src/dimensional_traits_3.h`.

Upon construction, the rhomboid tiling is computed from the provided input
points, up to the specified highest order. It provides access methods
for various combinatorial structures:

- The rhomboids themselves: `get_rhomboids`
- The order-k Delaunay mosaic: `get_slice_mosaic`
- The degree-k Delaunay mosaic: `get_halfint_slice_mosaic`
- Delaunay _slabs_: `get_slab_mosaic`

The order-k Delaunay mosaic is obtained as the intersection of the rhomboid
tiling with a hyperplane at depth k; degree-k Delaunay mosaic as the
intersection with a hyperplane at depth k+0.5; a Delaunay slab is the
intersection of the space between the hyperplanes at depth k and k+1 and
the rhomboid tiling.

The rhomboid tiling is equipped with a radius function on its rhomboids,
which allows for the computation of a discrete version of the 2-parameter
filtration (radius r and depth k) of multi-covers of Euclidean balls.
The related filtration on the order-k Delaunay mosaics is a discrete version
of the filtration by parameter r of k-covers of Euclidean balls for a fixed k.
Access methods for these filtration are:

- Filtration of the rhomboid tiling: `get_rhomboid_filtration`
- Filtration of the order-k Delaunay mosaic: `get_slice_filtration` and `get_delaunay_filtration`
- 2-parameter filtration of sliced rhomboids: `get_bifiltration`
- 2-parameter filtration of unsliced rhomboids: `get_rhomboid_bifiltration`

See the documentation in `src/rhomboid_tiling.h` for details.


## Commandline tool

The commandline tool accepts an input filename, output filename, dimension
parameter, order up to which to compute the mosaics, and an optional parameter
indicating the desired output data.

Example usage:
```
./main example_data/example_input.txt example_data/example_output.txt 3 4 slices
```
Input and output are text files. The input is one point per line,
each with d coordinates separated by spaces, where d is the provided dimension
parameter (currently only 2 and 3 dimensional data is supported).
Options for the output data are:

- bifi: [default] Boundary matrix of the 2-parameter discrete sliced rhomboid bifiltration.
- cbifi: A combinatorial representation of bifi (a cell is represented by its vertices rather than its ID).
- ubifi: boundary matrix of (unsliced) rhomboid bifiltration.
- rhomboids: The rhomboids of the rhomboid tiling (top dimensional cells only).
- slices: All order-k Delaunay slices (top dimensional cells only).
- halfint: All degree-k Delaunay half-integer slices (top dimensional cells only).
- slabs: All Delaunay slabs, i.e. cells between order-k and order-(k+1) Delaunay slices (top dimensional cells only).
- frhomboids: The filtration on the rhomboid tiling.
- fslices: Filtrations on the order-k Delaunay slices for each k.
- fslabs: Filtrations on the order-k Delaunay slices for each k and intermediate Delaunay slabs.
- bslices: Boundary matrices for order-k Delaunay slice filtration for each k.
- firep: free implicit representation of bifi.
- ufirep: free implicit representation of ubifi.

The output format depends on the output data. For bifi and ubifi, each cell is
assigned an ID, and in the output file, each cell is line of space separated
values, consisting of its ID, its dimension, its order *k*, its radius value
*r*, and the IDs of its boundary cells.
For the other output formats, the output meaning and data type
are documented in `src/rhomboid_tiling.h`. The formatting of each output data
type is defined in the respective class definition, with the exception of
CCell, which is a `std::vector` of `std::vector` of `int`s and is output in 
the list format used in python, i.e. comma separated values wrapped in square
brackets.
