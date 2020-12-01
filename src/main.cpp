/*
 * Copyright (c) 2019-2020 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "dimensional_traits_3.h"
#include "dimensional_traits_2.h"
#include "rhomboid_tiling.h"
#include "rhomboid.h"
#include "utils.h"

#include "cell_with_radius.h"
#include "rhomboid_with_radius.h"
#include "bifiltration_cell.h"
#include "combinatorial_bifiltration_cell.h"
#include "firep.h"

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Exact_predicates_exact_constructions_kernel                    K;
typedef K::FT                                                               FT;
typedef DimensionalTraits_3<K>                                             Dt3;
typedef DimensionalTraits_2<K>                                             Dt2;
typedef CombinatorialBifiltrationCell<FT>                                  CBC;


// Read points from a file and return a vector of points.
// File should consist only of points of dimension Dt::dimension;
// Each point is encoded as its space-separated coordinates.
template<class Dt>
std::vector<typename Dt::Point> read_points(std::istream &stream) {
    std::vector<typename Dt::Point> points;
    for (std::string s; std::getline(stream, s);) {
        if (s.size() == 0 || s[0] == '#')
            continue;

        std::stringstream ss(s);
        std::vector<double> p;
        double v;

        while (!ss.eof()) {
            ss >> v;
            if (ss.fail()) break;
            p.push_back(v);
        }
        if (p.size() == Dt::dimension) {
            points.push_back(Dt::make_point(p));
        } else if (!p.size() == 0) {
            // Line was not empty, but also not valid point.
            // --> Indicate invalid input
            return std::vector<typename Dt::Point>();
        }
    }
    return points;
}


void print_usage() {
    std::cout << "Usage: ./orderk infile outfile dimension order [type]" << std::endl;
    std::cout << "infile: Text file with 3 space separated coordinates per line." << std::endl;
    std::cout << "outfile: Output filename." << std::endl;
    std::cout << "dimension: Dimension of the input point set (2 or 3)." << std::endl;
    std::cout << "order: Order k up to which to compute requested structure." << std::endl;
    std::cout << "type: One of the following:" << std::endl;
    std::cout << "  bifi: [default] boundary matrix of (sliced) rhomboid bifiltration." << std::endl;
    std::cout << "  cbifi: combinatorial representation of bifi." << std::endl;
    std::cout << "  ubifi: boundary matrix of unsliced rhomboid bifiltration." << std::endl;
    std::cout << "  rhomboids: rhomboid tiling (top dimensional cells)." << std::endl;
    std::cout << "  slices: all delaunay slices (top dimensional cells)." << std::endl;
    std::cout << "  halfint: all delaunay half-integer slices (top dimensional cells)." << std::endl;
    std::cout << "  slabs: all delaunay slabs (top dimensional cells)." << std::endl;
    std::cout << "  frhomboids: rhomboid filtration." << std::endl;
    std::cout << "  fslices: slice filtrations." << std::endl;
    std::cout << "  fslabs: slice and slab filtrations." << std::endl;
    std::cout << "  bslices: boundary matrix for each slice filtration." << std::endl;
    std::cout << "  firep: free implicit representation of bifi." << std::endl;
    std::cout << "  ufirep: free implicit representation of ubifi." << std::endl;
}


template<class Dt>
void process_request(std::ifstream& pfile, std::ofstream& ofile, int max_order, std::string otype, int repr_dimension) {

    auto points = read_points<Dt>(pfile);
    pfile.close();
    if (points.size() == 0) {
        std::cout << "invalid input format" << std::endl;
        exit(0); 
    }

    std::cout << "Points loaded." << std::endl;

    if (max_order > points.size()) {
        max_order = points.size() - 1;
    }

    auto rt = RhomboidTiling<Dt>(points, max_order);
    std::cout << "Writing output." << std::endl;
    if (otype == "cbifi") {
        auto bf = rt.get_bifiltration();
        auto cm = rt.get_bifiltration_id_map();
        std::vector<CBC> cbcs;
        for (const auto& bc : bf) {
            cbcs.push_back(CBC(bc, cm));
        }
        std::sort(cbcs.begin(), cbcs.end());
        for (const CBC& c : cbcs) {
            ofile << c << std::endl;
        }
    } else if (otype == "rhomboids") {
        auto rhos = rt.get_rhomboids();
        for (const auto& rho : rhos) {
            ofile << rho << std::endl;
        }
    } else if (otype == "slices") {
        for (int order = 1; order <= max_order; ++order) {
            auto cells = rt.get_slice_mosaic(order);
            std::sort(cells.begin(), cells.end());
            ofile << "Slice " << order << ":" << std::endl;
            for (const auto& cell : cells) {
                ofile << cell << std::endl;
            }
        }
    } else if (otype == "halfint") {
        for (int order = 1; order < max_order; ++order) {
            auto cells = rt.get_halfint_slice_mosaic(order);
            std::sort(cells.begin(), cells.end());
            ofile << "Slice " << order+0.5 << ":" << std::endl;
            for (const auto& cell : cells) {
                ofile << cell << std::endl;
            }
        }
    } else if (otype == "slabs") {
        for (int order = 0; order < max_order; ++order) {
            auto cells = rt.get_slab_mosaic(order);
            std::sort(cells.begin(), cells.end());
            ofile << "Slab between " << order << " and " << order+1 << ":" << std::endl;
            for (const auto& cell : cells) {
                ofile << cell << std::endl;
            }
        }
    } else if (otype == "frhomboids") {
        auto rhos = rt.get_rhomboid_filtration();
        for (const auto& rho : rhos) {
            ofile << rho << std::endl;
        }
    } else if (otype == "fslices") {
        for (int order = 1; order <= max_order; ++order) {
            ofile << "Slice " << order << ":" << std::endl;
            auto cells = rt.get_slice_filtration(order);
            for (const auto& cell : cells) {
                ofile << cell << std::endl;
            }
        }
    } else if (otype == "fslabs") {
        for (int order = 1; order < max_order; ++order) {
            ofile << "Slice " << order << ":" << std::endl;
            auto cells = rt.get_slice_filtration(order);
            for (const auto& cell : cells) {
                ofile << cell << std::endl;
            }
            ofile << "Slab between " << order << " and " << order+1 <<  ":" << std::endl;
            cells = rt.get_slab_filtration(order);
            for (const auto& cell : cells) {
                ofile << cell << std::endl;
            }
        }
        ofile << "Slice " << max_order << ":" << std::endl;
        auto cells = rt.get_slice_filtration(max_order);
        for (const auto& cell : cells) {
            ofile << cell << std::endl;
        }
    } else if (otype == "bslices") {
        for (int order = 1; order < max_order; ++order) {
            ofile << "Slice " << order << ":" << std::endl;
            auto bf = rt.get_delaunay_filtration(order);
            std::sort(bf.begin(), bf.end());
            for (const auto& c : bf) {
                ofile << c << std::endl;
            }
        }
    } else if (otype == "ubifi") {
        auto bf = rt.get_unsliced_bifiltration();
        std::sort(bf.begin(), bf.end());
        for (const auto& c : bf) {
            ofile << c << std::endl;
        }
    } else if (otype == "ufirep") {
        auto bf = rt.get_unsliced_bifiltration();
        bifiltration_to_firep<FT>(bf, repr_dimension, ofile);
    } else if (otype == "firep") {
        auto bf = rt.get_bifiltration();
        bifiltration_to_firep<FT>(bf, repr_dimension, ofile);
    } else if (otype == "bifi") { // "bifi" is the default
        auto bf = rt.get_bifiltration();
        std::sort(bf.begin(), bf.end());
        for (const auto& c : bf) {
            ofile << c << std::endl;
        }
    } else {
        std::cout << "Invalid output type specified, aborting." << std::endl;
    }
}


int main(int argc, char** argv)
{
    if (argc < 5) {
        print_usage();
        exit(0);
    }

    std::string infile = std::string(argv[1]);
    std::string outfile = std::string(argv[2]);
    int dimension = std::atoi(argv[3]);
    int max_order = std::atoi(argv[4]);
    std::string otype = "bifi";
    if (argc > 5) {
        otype = std::string(argv[5]);
    }
    // for firep only
    int repr_dimension = 1;
    if (argc > 6) {
        repr_dimension = atoi(argv[6]);
    }

    std::ifstream pfile(infile.c_str());
    if (!pfile) {
        std::cout << "infile not found" << std::endl;
        exit(0);
    }

    std::ofstream ofile(outfile.c_str());
    if (!ofile) {
        std::cout << "outfile not found" << std::endl;
        exit(0);
    }

    if (max_order <= 0) { // we get 0 if the atoi conversion fails.
        std::cout << "invalid order" << std::endl;
        exit(0);
    }

    if (dimension == 2) {
        process_request<Dt2>(pfile, ofile, max_order, otype, repr_dimension);
    } else if (dimension == 3) {
        process_request<Dt3>(pfile, ofile, max_order, otype, repr_dimension);
    } else {
        std::cout << "invalid dimension" << std::endl;
        exit(0);        
    }

    ofile.close();
    std::cout << "DONE." << std::endl;
}
