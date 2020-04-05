/*
 * Copyright (c) 2019 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#include "dimensional_traits_3.h"
#include "rhomboid_tiling.h"
#include "rhomboid.h"
#include "utils.h"

#include "cell_with_radius.h"
#include "rhomboid_with_radius.h"
#include "bifiltration_cell.h"
#include "combinatorial_bifiltration_cell.h"

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

using namespace Catch;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "testutils.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel                    K;
typedef K::FT                                                               FT;
typedef DimensionalTraits_3<K>                                             Dt3;
typedef RhomboidTiling<Dt3>                                   RhomboidTiling_3;
typedef RhomboidTiling_3::Point                                          Point;
typedef CombinatorialBifiltrationCell<FT>                                  CBC;


TEST_CASE("Standard toy example for rhomboids", "[rhombids][toy]") {

    Point p0(0, 0, 0);
    Point p1(0, 4, 4);
    Point p2(4, 4, 0);
    Point p3(4, 0, 4);
    Point p4(-10, 2, 2);

    std::vector<Point> points = {p0, p1, p2, p3, p4};

    std::vector<Rhomboid> rhos_expected = {
        Rhomboid({}, {0, 1, 2, 4}),
        Rhomboid({}, {0, 1, 2, 3}),
        Rhomboid({}, {0, 1, 3, 4}),
        Rhomboid({0}, {1, 2, 3, 4}),
        Rhomboid({1}, {0, 2, 3, 4}),
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 5);
    auto rhos = rhomboidtiling.get_rhomboids();
    CHECK_THAT(rhos, UnorderedEquals(rhos_expected));
}


TEST_CASE("Standard toy example for orderk", "[orderk][toy]") {

    Point p0(0, 0, 0);
    Point p1(0, 4, 4);
    Point p2(4, 4, 0);
    Point p3(4, 0, 4);
    Point p4(-10, 2, 2);

    std::vector<Point> points = {p0, p1, p2, p3, p4};

    std::vector<std::vector<CVertex>> o1del_expected = {
        {{0}, {1}, {2}, {3}},
        {{0}, {1}, {2}, {4}},
        {{0}, {1}, {3}, {4}}
    };

    std::vector<std::vector<CVertex>> o2del_expected = {
        {{0, 1}, {0, 2}, {0, 3}, {0, 4}},
        {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}},
        {{0, 1}, {0, 2}, {0, 4}, {1, 2}, {1, 4}, {2, 4}},
        {{0, 1}, {0, 3}, {0, 4}, {1, 3}, {1, 4}, {3, 4}},
        {{0, 1}, {1, 2}, {1, 3}, {1, 4}}
    };

    std::vector<std::vector<CVertex>> o3del_expected = {
        {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {0, 2, 3}, {0, 2, 4}, {0, 3, 4}},
        {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {1, 2, 3}, {1, 2, 4}, {1, 3, 4}},
        {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}},
        {{0, 1, 2}, {0, 1, 4}, {0, 2, 4}, {1, 2, 4}},
        {{0, 1, 3}, {0, 1, 4}, {0, 3, 4}, {1, 3, 4}}
    };

    std::vector<std::vector<CVertex>> o4del_expected = {
        {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}},
        {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {1, 2, 3, 4}}
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 4);

    auto o1del = rhomboidtiling.get_slice_mosaic(1);
    CHECK_THAT(o1del, UnorderedEquals(o1del_expected));

    auto o2del = rhomboidtiling.get_slice_mosaic(2);
    CHECK_THAT(o2del, UnorderedEquals(o2del_expected));

    auto o3del = rhomboidtiling.get_slice_mosaic(3);
    CHECK_THAT(o3del, UnorderedEquals(o3del_expected));

    auto o4del = rhomboidtiling.get_slice_mosaic(4);
    CHECK_THAT(o4del, UnorderedEquals(o4del_expected));
}


TEST_CASE("Minimal example with a single cell", "[orderk][minimal]") {

    Point p0(0, 0, 0);
    Point p1(0, 4, 4);
    Point p2(4, 4, 0);
    Point p3(4, 0, 4);

    std::vector<Point> points = {p0, p1, p2, p3};

    std::vector<std::vector<CVertex>> o1del_expected = {
        {{0,}, {1,}, {2,}, {3,}}
    };

    std::vector<std::vector<CVertex>> o2del_expected = {
        {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}
    };

    std::vector<std::vector<CVertex>> o3del_expected = {
        {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 3);

    auto o1del = rhomboidtiling.get_slice_mosaic(1);
    CHECK_THAT(o1del, UnorderedEquals(o1del_expected));

    auto o2del = rhomboidtiling.get_slice_mosaic(2);
    CHECK_THAT(o2del, UnorderedEquals(o2del_expected));

    auto o3del = rhomboidtiling.get_slice_mosaic(3);
    CHECK_THAT(o3del, UnorderedEquals(o3del_expected));
}


TEST_CASE("Non-convex cluster", "[orderk][cluster]") {

    Point p0(0,0,3);
    Point p1(0,-1.5,-1.8);
    Point p2(-0.07,3.67,-2.03);
    Point p3(-2.37,3.08,2.49);
    Point p4(2.32,4.37,0.4);
    Point p5(0,-1.5,0);

    std::vector<Point> points = {p0, p1, p2, p3, p4, p5};

    // The cells
    // {{0, 5}, {1, 5}, {2, 5}, {3, 5}} and
    // {{0, 5}, {1, 5}, {2, 5}, {4, 5}}
    // form a non-convex cluster.
    std::vector<std::vector<CVertex>> o2del_expected = {
        {{0, 2}, {0, 3}, {0, 4}, {0, 5}},
        {{0, 2}, {0, 3}, {0, 4}, {2, 3}, {2, 4}, {3, 4}},
        {{0, 2}, {0, 3}, {0, 5}, {2, 3}, {2, 5}, {3, 5}},
        {{0, 2}, {0, 4}, {0, 5}, {2, 4}, {2, 5}, {4, 5}},
        {{0, 2}, {2, 3}, {2, 4}, {2, 5}},
        {{0, 4}, {1, 4}, {2, 4}, {4, 5}},
        {{0, 5}, {1, 5}, {2, 5}, {3, 5}},
        {{0, 5}, {1, 5}, {2, 5}, {4, 5}},
        {{1, 2}, {1, 3}, {1, 5}, {2, 3}, {2, 5}, {3, 5}},
        {{1, 2}, {1, 4}, {1, 5}, {2, 4}, {2, 5}, {4, 5}},
        {{1, 2}, {2, 3}, {2, 4}, {2, 5}}
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 2);

    auto o2del = rhomboidtiling.get_slice_mosaic(2);
    CHECK_THAT(o2del, UnorderedEquals(o2del_expected));
}


TEST_CASE("Radius for noncritical rhomboid", "[orderkfilt][minimal][noncritical]") {

    Point p0(0, 0, 0);
    Point p1(0,-1, 2);
    Point p2(0, 2,-1);
    Point p3(1, 0, 0);

    std::vector<Point> points = {p0, p1, p2, p3};

    std::vector<CellWithRadius<FT>> o1del_expected = {
        CellWithRadius<FT>({{0}}, 0),
        CellWithRadius<FT>({{1}}, 0),
        CellWithRadius<FT>({{2}}, 0),
        CellWithRadius<FT>({{3}}, 0),
        CellWithRadius<FT>({{0}, {1}}, 1.25),
        CellWithRadius<FT>({{0}, {2}}, 1.25),
        CellWithRadius<FT>({{0}, {3}}, 0.25),
        CellWithRadius<FT>({{1}, {2}}, 12.5),
        CellWithRadius<FT>({{1}, {3}}, 1.5),
        CellWithRadius<FT>({{2}, {3}}, 1.5),
        CellWithRadius<FT>({{0}, {1}, {2}}, 12.5),
        CellWithRadius<FT>({{0}, {1}, {3}}, 1.5),
        CellWithRadius<FT>({{0}, {2}, {3}}, 1.5),
        CellWithRadius<FT>({{1}, {2}, {3}}, 12.75),
        CellWithRadius<FT>({{0}, {1}, {2}, {3}}, 12.75)
    };

    std::vector<CellWithRadius<FT>> o2del_expected = {
        CellWithRadius<FT>({{0, 1}}, 1.25),
        CellWithRadius<FT>({{0, 2}}, 1.25),
        CellWithRadius<FT>({{0, 3}}, 0.25),
        CellWithRadius<FT>({{1, 2}}, 12.5),
        CellWithRadius<FT>({{1, 3}}, 1.5),
        CellWithRadius<FT>({{2, 3}}, 1.5),
        CellWithRadius<FT>({{0, 1}, {0, 2}}, 6.0),
        CellWithRadius<FT>({{0, 1}, {0, 3}}, 1.5),
        CellWithRadius<FT>({{0, 1}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{0, 1}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{0, 2}, {0, 3}}, 1.5),
        CellWithRadius<FT>({{0, 2}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{0, 2}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{0, 3}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{0, 3}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{1, 2}, {1, 3}}, 12.75),
        CellWithRadius<FT>({{1, 2}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{1, 3}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {0, 3}}, 6.0),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{0, 1}, {0, 3}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{0, 1}, {1, 2}, {1, 3}}, 12.75),
        CellWithRadius<FT>({{0, 2}, {0, 3}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{0, 2}, {1, 2}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 3}, {1, 3}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{1, 2}, {1, 3}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 12.75)
    };

    std::vector<CellWithRadius<FT>> o3del_expected = {
        CellWithRadius<FT>({{0, 1, 2}}, 6.0),
        CellWithRadius<FT>({{0, 1, 3}}, 1.5),
        CellWithRadius<FT>({{0, 2, 3}}, 1.5),
        CellWithRadius<FT>({{1, 2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}}, 6.0),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}}, 6.0),
        CellWithRadius<FT>({{0, 1, 2}, {1, 2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1, 3}, {0, 2, 3}}, 4.5),
        CellWithRadius<FT>({{0, 1, 3}, {1, 2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 2, 3}, {1, 2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {0, 2, 3}}, 6.0),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {1, 2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}, {1, 2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1, 3}, {0, 2, 3}, {1, 2, 3}}, 12.75),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}, 12.75)
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 4);
    auto o1del = rhomboidtiling.get_slice_filtration(1);
    auto o2del = rhomboidtiling.get_slice_filtration(2);
    auto o3del = rhomboidtiling.get_slice_filtration(3);
    // // TODO: Figure out why these rounding errors occur.
    // for (int i = 0; i < o1del.size(); ++i) {
    //     CHECK(o1del[i].r == o1del_expected[i].r);
    // }
    CHECK_THAT(o1del, UnorderedEquals(o1del_expected));
    CHECK_THAT(o2del, UnorderedEquals(o2del_expected));
    CHECK_THAT(o3del, UnorderedEquals(o3del_expected));
}


// TODO: make values exact
TEST_CASE("Radius for toy example", "[orderkfilt][toy]") {

    Point p0(0, 0, 0);
    Point p1(0, 4, 4);
    Point p2(4, 4, 0);
    Point p3(4, 0, 4);
    Point p4(-10, 2, 2);

    std::vector<Point> points = {p0, p1, p2, p3, p4};

    std::vector<CellWithRadius<FT>> o1del_expected = {
        CellWithRadius<FT>({{0}}, 0),
        CellWithRadius<FT>({{1}}, 0),
        CellWithRadius<FT>({{2}}, 0),
        CellWithRadius<FT>({{3}}, 0),
        CellWithRadius<FT>({{4}}, 0),
        CellWithRadius<FT>({{0}, {1}}, 8),
        CellWithRadius<FT>({{0}, {2}}, 8),
        CellWithRadius<FT>({{0}, {3}}, 8),
        CellWithRadius<FT>({{0}, {4}}, 27),
        CellWithRadius<FT>({{1}, {2}}, 8),
        CellWithRadius<FT>({{1}, {3}}, 8),
        CellWithRadius<FT>({{1}, {4}}, 27),
        CellWithRadius<FT>({{2}, {3}}, 8),
        CellWithRadius<FT>({{2}, {4}}, 116.28),
        CellWithRadius<FT>({{3}, {4}}, 116.28),
        CellWithRadius<FT>({{0}, {1}, {2}}, 10.6666666667),
        CellWithRadius<FT>({{0}, {1}, {3}}, 10.6666666667),
        CellWithRadius<FT>({{0}, {1}, {4}}, 29.16),
        CellWithRadius<FT>({{0}, {2}, {3}}, 10.6666666667),
        CellWithRadius<FT>({{0}, {2}, {4}}, 116.28),
        CellWithRadius<FT>({{0}, {3}, {4}}, 116.28),
        CellWithRadius<FT>({{1}, {2}, {3}}, 10.6666666667),
        CellWithRadius<FT>({{1}, {2}, {4}}, 116.28),
        CellWithRadius<FT>({{1}, {3}, {4}}, 116.28),
        CellWithRadius<FT>({{0}, {1}, {2}, {3}}, 12),
        CellWithRadius<FT>({{0}, {1}, {2}, {4}}, 116.28),
        CellWithRadius<FT>({{0}, {1}, {3}, {4}}, 116.28)
    };

    std::vector<CellWithRadius<FT>> o2del_expected = {
        CellWithRadius<FT>({{0, 1}}, 8),
        CellWithRadius<FT>({{0, 2}}, 8),
        CellWithRadius<FT>({{0, 3}}, 8),
        CellWithRadius<FT>({{0, 4}}, 27),
        CellWithRadius<FT>({{1, 2}}, 8),
        CellWithRadius<FT>({{1, 3}}, 8),
        CellWithRadius<FT>({{1, 4}}, 27),
        CellWithRadius<FT>({{2, 3}}, 8),
        CellWithRadius<FT>({{2, 4}}, 116.28),
        CellWithRadius<FT>({{3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1}, {0, 2}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1}, {0, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1}, {0, 4}}, 29.16),
        CellWithRadius<FT>({{0, 1}, {1, 2}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1}, {1, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1}, {1, 4}}, 29.16),
        CellWithRadius<FT>({{0, 2}, {0, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 2}, {0, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 2}, {1, 2}}, 10.6666666667),
        CellWithRadius<FT>({{0, 2}, {2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 2}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 3}, {0, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 3}, {1, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 3}, {2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 3}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 4}, {1, 4}}, 29.16),
        CellWithRadius<FT>({{0, 4}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 4}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{1, 2}, {1, 3}}, 10.6666666667),
        CellWithRadius<FT>({{1, 2}, {1, 4}}, 72.4736842105),
        CellWithRadius<FT>({{1, 2}, {2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{1, 2}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{1, 3}, {1, 4}}, 72.4736842105),
        CellWithRadius<FT>({{1, 3}, {2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{1, 3}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{1, 4}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{1, 4}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {0, 3}}, 12),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {0, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {1, 2}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1}, {0, 3}, {0, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1}, {0, 3}, {1, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1}, {0, 4}, {1, 4}}, 29.16),
        CellWithRadius<FT>({{0, 1}, {1, 2}, {1, 3}}, 12),
        CellWithRadius<FT>({{0, 1}, {1, 2}, {1, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1}, {1, 3}, {1, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 2}, {0, 3}, {0, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 2}, {0, 3}, {2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 2}, {0, 4}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 2}, {1, 2}, {2, 3}}, 12),
        CellWithRadius<FT>({{0, 2}, {1, 2}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 3}, {0, 4}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 3}, {1, 3}, {2, 3}}, 12),
        CellWithRadius<FT>({{0, 3}, {1, 3}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 4}, {1, 4}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 4}, {1, 4}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{1, 2}, {1, 3}, {1, 4}}, 97.5306122449),
        CellWithRadius<FT>({{1, 2}, {1, 3}, {2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{1, 2}, {1, 4}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{1, 3}, {1, 4}, {3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {0, 3}, {0, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1}, {1, 2}, {1, 3}, {1, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 12),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {0, 4}, {1, 2}, {1, 4}, {2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1}, {0, 3}, {0, 4}, {1, 3}, {1, 4}, {3, 4}}, 116.28)
    };

    std::vector<CellWithRadius<FT>> o3del_expected = {
        CellWithRadius<FT>({{0, 1, 2}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 1, 4}}, 29.16),
        CellWithRadius<FT>({{0, 2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{0, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{1, 2, 3}}, 10.6666666667),
        CellWithRadius<FT>({{1, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{1, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 4}}, 51),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 2}, {1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2}, {1, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 3}, {0, 1, 4}}, 51),
        CellWithRadius<FT>({{0, 1, 3}, {0, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 3}, {0, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 3}, {1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 3}, {1, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 4}, {0, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 4}, {0, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 4}, {1, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 4}, {1, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 2, 3}, {0, 2, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 2, 3}, {0, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 2, 3}, {1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 2, 4}, {0, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 2, 4}, {1, 2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 3, 4}, {1, 3, 4}}, 116.28),
        CellWithRadius<FT>({{1, 2, 3}, {1, 2, 4}}, 97.5306122449),
        CellWithRadius<FT>({{1, 2, 3}, {1, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{1, 2, 4}, {1, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {0, 1, 4}}, 53.0816326531),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {0, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 4}, {0, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 4}, {1, 2, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}, {0, 2, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}, {1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 4}, {1, 2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1, 2}, {1, 2, 3}, {1, 2, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 3}, {0, 1, 4}, {0, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 3}, {0, 1, 4}, {1, 3, 4}}, 72.4736842105),
        CellWithRadius<FT>({{0, 1, 3}, {0, 2, 3}, {0, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 3}, {0, 2, 3}, {1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 3}, {0, 3, 4}, {1, 3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1, 3}, {1, 2, 3}, {1, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 4}, {0, 2, 4}, {0, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 4}, {0, 2, 4}, {1, 2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1, 4}, {0, 3, 4}, {1, 3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1, 4}, {1, 2, 4}, {1, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 2, 3}, {0, 2, 4}, {0, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{1, 2, 3}, {1, 2, 4}, {1, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 4}, {0, 2, 4}, {1, 2, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1, 3}, {0, 1, 4}, {0, 3, 4}, {1, 3, 4}}, 116.28),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {0, 2, 3}, {0, 2, 4}, {0, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {1, 2, 3}, {1, 2, 4}, {1, 3, 4}}, 97.5306122449)
    };

    std::vector<CellWithRadius<FT>> o4del_expected = {
        CellWithRadius<FT>({{0, 1, 2, 3}}, 12),
        CellWithRadius<FT>({{0, 1, 2, 4}}, 51),
        CellWithRadius<FT>({{0, 1, 3, 4}}, 51),
        CellWithRadius<FT>({{0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{1, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}}, 53.0816326531),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 3, 4}}, 53.0816326531),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {1, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 4}, {0, 1, 3, 4}}, 53.0816326531),
        CellWithRadius<FT>({{0, 1, 2, 4}, {0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 4}, {1, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 3, 4}, {0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 3, 4}, {1, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}}, 53.0816326531),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}, {1, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 3, 4}, {0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 3, 4}, {1, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 4}, {0, 1, 3, 4}, {1, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}}, 97.5306122449),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {1, 2, 3, 4}}, 97.5306122449)
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 5);
    auto o1del = rhomboidtiling.get_slice_filtration(1);
    auto o2del = rhomboidtiling.get_slice_filtration(2);
    auto o3del = rhomboidtiling.get_slice_filtration(3);
    auto o4del = rhomboidtiling.get_slice_filtration(4);
    CHECK_THAT(o1del, UnorderedEquals(o1del_expected));
    CHECK_THAT(o2del, UnorderedEquals(o2del_expected));
    CHECK_THAT(o3del, UnorderedEquals(o3del_expected));
    CHECK_THAT(o4del, UnorderedEquals(o4del_expected));
}



TEST_CASE("Slabs of non-critical rhomboid", "[orderkslab][minimal]") {

    Point p0(0, 0, 0);
    Point p1(0,-1, 2);
    Point p2(0, 2,-1);
    Point p3(1, 0, 0);

    std::vector<Point> points = {p0, p1, p2, p3};

    std::vector<CellWithRadius<FT>> o01del_expected = {
        CellWithRadius<FT>({{}, {0}}, 0),
        CellWithRadius<FT>({{}, {1}}, 0),
        CellWithRadius<FT>({{}, {2}}, 0),
        CellWithRadius<FT>({{}, {3}}, 0),
        CellWithRadius<FT>({{}, {0}, {1}}, 1.25),
        CellWithRadius<FT>({{}, {0}, {2}}, 1.25),
        CellWithRadius<FT>({{}, {0}, {3}}, 0.25),
        CellWithRadius<FT>({{}, {1}, {2}}, 12.5),
        CellWithRadius<FT>({{}, {1}, {3}}, 1.5),
        CellWithRadius<FT>({{}, {2}, {3}}, 1.5),
        CellWithRadius<FT>({{}, {0}, {1}, {2}}, 12.5),
        CellWithRadius<FT>({{}, {0}, {1}, {3}}, 1.5),
        CellWithRadius<FT>({{}, {0}, {2}, {3}}, 1.5),
        CellWithRadius<FT>({{}, {1}, {2}, {3}}, 12.75),
        CellWithRadius<FT>({{}, {0}, {1}, {2}, {3}}, 12.75),
    };

    std::vector<CellWithRadius<FT>> o12del_expected = {
        CellWithRadius<FT>({{0}, {0, 1}}, 1.25),
        CellWithRadius<FT>({{0}, {0, 2}}, 1.25),
        CellWithRadius<FT>({{0}, {0, 3}}, 0.25),
        CellWithRadius<FT>({{1}, {0, 1}}, 1.25),
        CellWithRadius<FT>({{1}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{1}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{2}, {0, 2}}, 1.25),
        CellWithRadius<FT>({{2}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{2}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{3}, {0, 3}}, 0.25),
        CellWithRadius<FT>({{3}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{3}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{0}, {1}, {0, 1}}, 1.25),
        CellWithRadius<FT>({{0}, {2}, {0, 2}}, 1.25),
        CellWithRadius<FT>({{0}, {3}, {0, 3}}, 0.25),
        CellWithRadius<FT>({{1}, {2}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{1}, {3}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{2}, {3}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{0}, {0, 1}, {0, 2}}, 6),
        CellWithRadius<FT>({{0}, {0, 1}, {0, 3}}, 1.5),
        CellWithRadius<FT>({{0}, {0, 2}, {0, 3}}, 1.5),
        CellWithRadius<FT>({{1}, {0, 1}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{1}, {0, 1}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{1}, {1, 2}, {1, 3}}, 12.75),
        CellWithRadius<FT>({{2}, {0, 2}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{2}, {0, 2}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{2}, {1, 2}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{3}, {0, 3}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{3}, {0, 3}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{3}, {1, 3}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{0}, {0, 1}, {0, 2}, {0, 3}}, 6),
        CellWithRadius<FT>({{1}, {0, 1}, {1, 2}, {1, 3}}, 12.75),
        CellWithRadius<FT>({{2}, {0, 2}, {1, 2}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{3}, {0, 3}, {1, 3}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{0}, {1}, {2}, {0, 1}, {0, 2}, {1, 2}}, 12.5),
        CellWithRadius<FT>({{0}, {1}, {3}, {0, 1}, {0, 3}, {1, 3}}, 1.5),
        CellWithRadius<FT>({{0}, {2}, {3}, {0, 2}, {0, 3}, {2, 3}}, 1.5),
        CellWithRadius<FT>({{1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}}, 12.75),
        CellWithRadius<FT>({{0}, {1}, {2}, {3}, {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 12.75),
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 2);
    auto o01del = rhomboidtiling.get_slab_filtration(0);
    auto o12del = rhomboidtiling.get_slab_filtration(1);
    CHECK_THAT(o01del, UnorderedEquals(o01del_expected));
    CHECK_THAT(o12del, UnorderedEquals(o12del_expected));
}


TEST_CASE("Bifiltration of non-critical rhomboid", "[bifi][minimal]") {
    // We cannot directly compare the bifiltrations, because the IDs
    // assigned to cells are not unique. Thus we expand the IDs into
    // combinatorial cells, sort the output and compare.

    Point p0(0, 0, 0);
    Point p1(0,-1, 2);
    Point p2(0, 2,-1);
    Point p3(1, 0, 0);

    std::vector<Point> points = {p0, p1, p2, p3};

    const std::vector<CBC> cbifi_expected = {
        CBC(0, 1, {{0}}, 0, {}),
        CBC(0, 1, {{1}}, 0, {}),
        CBC(0, 1, {{2}}, 0, {}),
        CBC(0, 1, {{3}}, 0, {}),
        CBC(0, 2, {{0, 1}}, 1.25, {}),
        CBC(0, 2, {{0, 2}}, 1.25, {}),
        CBC(0, 2, {{0, 3}}, 0.25, {}),
        CBC(0, 2, {{1, 2}}, 12.5, {}),
        CBC(0, 2, {{1, 3}}, 1.5, {}),
        CBC(0, 2, {{2, 3}}, 1.5, {}),
        CBC(1, 1, {{0}, {0, 1}}, 1.25, {{{0}}, {{0, 1}}}),
        CBC(1, 1, {{0}, {0, 2}}, 1.25, {{{0}}, {{0, 2}}}),
        CBC(1, 1, {{0}, {0, 3}}, 0.25, {{{0}}, {{0, 3}}}),
        CBC(1, 1, {{0}, {1}}, 1.25, {{{0}}, {{1}}}),
        CBC(1, 1, {{0}, {2}}, 1.25, {{{0}}, {{2}}}),
        CBC(1, 1, {{0}, {3}}, 0.25, {{{0}}, {{3}}}),
        CBC(1, 1, {{0, 1}, {1}}, 1.25, {{{1}}, {{0, 1}}}),
        CBC(1, 1, {{0, 2}, {2}}, 1.25, {{{2}}, {{0, 2}}}),
        CBC(1, 1, {{0, 3}, {3}}, 0.25, {{{3}}, {{0, 3}}}),
        CBC(1, 1, {{1}, {1, 2}}, 12.5, {{{1}}, {{1, 2}}}),
        CBC(1, 1, {{1}, {1, 3}}, 1.5, {{{1}}, {{1, 3}}}),
        CBC(1, 1, {{1}, {2}}, 12.5, {{{1}}, {{2}}}),
        CBC(1, 1, {{1}, {3}}, 1.5, {{{1}}, {{3}}}),
        CBC(1, 1, {{1, 2}, {2}}, 12.5, {{{2}}, {{1, 2}}}),
        CBC(1, 1, {{1, 3}, {3}}, 1.5, {{{3}}, {{1, 3}}}),
        CBC(1, 1, {{2}, {2, 3}}, 1.5, {{{2}}, {{2, 3}}}),
        CBC(1, 1, {{2}, {3}}, 1.5, {{{2}}, {{3}}}),
        CBC(1, 1, {{2, 3}, {3}}, 1.5, {{{3}}, {{2, 3}}}),
        CBC(1, 2, {{0, 1}, {0, 2}}, 6, {{{0, 1}}, {{0, 2}}}),
        CBC(1, 2, {{0, 1}, {0, 3}}, 1.5, {{{0, 1}}, {{0, 3}}}),
        CBC(1, 2, {{0, 1}, {1, 2}}, 12.5, {{{0, 1}}, {{1, 2}}}),
        CBC(1, 2, {{0, 1}, {1, 3}}, 1.5, {{{0, 1}}, {{1, 3}}}),
        CBC(1, 2, {{0, 2}, {0, 3}}, 1.5, {{{0, 2}}, {{0, 3}}}),
        CBC(1, 2, {{0, 2}, {1, 2}}, 12.5, {{{0, 2}}, {{1, 2}}}),
        CBC(1, 2, {{0, 2}, {2, 3}}, 1.5, {{{0, 2}}, {{2, 3}}}),
        CBC(1, 2, {{0, 3}, {1, 3}}, 1.5, {{{0, 3}}, {{1, 3}}}),
        CBC(1, 2, {{0, 3}, {2, 3}}, 1.5, {{{0, 3}}, {{2, 3}}}),
        CBC(1, 2, {{1, 2}, {1, 3}}, 12.75, {{{1, 2}}, {{1, 3}}}),
        CBC(1, 2, {{1, 2}, {2, 3}}, 12.75, {{{1, 2}}, {{2, 3}}}),
        CBC(1, 2, {{1, 3}, {2, 3}}, 12.75, {{{1, 3}}, {{2, 3}}}),
        CBC(2, 1, {{0}, {0, 1}, {0, 2}}, 6, 
                  {{{0}, {0, 2}}, {{0}, {0, 1}}, {{0, 1}, {0, 2}}}),
        CBC(2, 1, {{0}, {0, 1}, {0, 3}}, 1.5, 
                  {{{0}, {0, 3}}, {{0}, {0, 1}}, {{0, 1}, {0, 3}}}),
        CBC(2, 1, {{0}, {0, 1}, {1}}, 1.25, 
                  {{{0}, {0, 1}}, {{0, 1}, {1}}, {{0}, {1}}}),
        CBC(2, 1, {{0}, {0, 2}, {0, 3}}, 1.5, 
                  {{{0}, {0, 3}}, {{0}, {0, 2}}, {{0, 2}, {0, 3}}}),
        CBC(2, 1, {{0}, {0, 2}, {2}}, 1.25, 
                  {{{0}, {0, 2}}, {{0, 2}, {2}}, {{0}, {2}}}),
        CBC(2, 1, {{0}, {0, 3}, {3}}, 0.25, 
                  {{{0}, {0, 3}}, {{0, 3}, {3}}, {{0}, {3}}}),
        CBC(2, 1, {{0}, {1}, {2}}, 12.5, 
                  {{{1}, {2}}, {{0}, {2}}, {{0}, {1}}}),
        CBC(2, 1, {{0}, {1}, {3}}, 1.5, 
                  {{{1}, {3}}, {{0}, {3}}, {{0}, {1}}}),
        CBC(2, 1, {{0}, {2}, {3}}, 1.5, 
                  {{{2}, {3}}, {{0}, {3}}, {{0}, {2}}}),
        CBC(2, 1, {{0, 1}, {1}, {1, 2}}, 12.5, 
                  {{{1}, {1, 2}}, {{0, 1}, {1}}, {{0, 1}, {1, 2}}}),
        CBC(2, 1, {{0, 1}, {1}, {1, 3}}, 1.5, 
                  {{{1}, {1, 3}}, {{0, 1}, {1}}, {{0, 1}, {1, 3}}}),
        CBC(2, 1, {{0, 2}, {1, 2}, {2}}, 12.5, 
                  {{{1, 2}, {2}}, {{0, 2}, {2}}, {{0, 2}, {1, 2}}}),
        CBC(2, 1, {{0, 2}, {2}, {2, 3}}, 1.5, 
                  {{{2}, {2, 3}}, {{0, 2}, {2}}, {{0, 2}, {2, 3}}}),
        CBC(2, 1, {{0, 3}, {1, 3}, {3}}, 1.5, 
                  {{{1, 3}, {3}}, {{0, 3}, {3}}, {{0, 3}, {1, 3}}}),
        CBC(2, 1, {{0, 3}, {2, 3}, {3}}, 1.5, 
                  {{{2, 3}, {3}}, {{0, 3}, {3}}, {{0, 3}, {2, 3}}}),
        CBC(2, 1, {{1}, {1, 2}, {1, 3}}, 12.75, 
                  {{{1}, {1, 3}}, {{1}, {1, 2}}, {{1, 2}, {1, 3}}}),
        CBC(2, 1, {{1}, {1, 2}, {2}}, 12.5, 
                  {{{1}, {1, 2}}, {{1, 2}, {2}}, {{1}, {2}}}),
        CBC(2, 1, {{1}, {1, 3}, {3}}, 1.5, 
                  {{{1}, {1, 3}}, {{1, 3}, {3}}, {{1}, {3}}}),
        CBC(2, 1, {{1}, {2}, {3}}, 12.75, 
                  {{{2}, {3}}, {{1}, {3}}, {{1}, {2}}}),
        CBC(2, 1, {{1, 2}, {2}, {2, 3}}, 12.75, 
                  {{{2}, {2, 3}}, {{1, 2}, {2}}, {{1, 2}, {2, 3}}}),
        CBC(2, 1, {{1, 3}, {2, 3}, {3}}, 12.75, 
                  {{{2, 3}, {3}}, {{1, 3}, {3}}, {{1, 3}, {2, 3}}}),
        CBC(2, 1, {{2}, {2, 3}, {3}}, 1.5, 
                  {{{2}, {2, 3}}, {{2, 3}, {3}}, {{2}, {3}}}),
        CBC(2, 2, {{0, 1}, {0, 2}, {0, 3}}, 6, 
                  {{{0, 2}, {0, 3}}, {{0, 1}, {0, 3}}, {{0, 1}, {0, 2}}}),
        CBC(2, 2, {{0, 1}, {0, 2}, {1, 2}}, 12.5, 
                  {{{0, 1}, {0, 2}}, {{0, 1}, {1, 2}}, {{0, 2}, {1, 2}}}),
        CBC(2, 2, {{0, 1}, {0, 3}, {1, 3}}, 1.5, 
                  {{{0, 1}, {0, 3}}, {{0, 1}, {1, 3}}, {{0, 3}, {1, 3}}}),
        CBC(2, 2, {{0, 1}, {1, 2}, {1, 3}}, 12.75, 
                  {{{1, 2}, {1, 3}}, {{0, 1}, {1, 3}}, {{0, 1}, {1, 2}}}),
        CBC(2, 2, {{0, 2}, {0, 3}, {2, 3}}, 1.5, 
                  {{{0, 2}, {0, 3}}, {{0, 2}, {2, 3}}, {{0, 3}, {2, 3}}}),
        CBC(2, 2, {{0, 2}, {1, 2}, {2, 3}}, 12.75, 
                  {{{1, 2}, {2, 3}}, {{0, 2}, {2, 3}}, {{0, 2}, {1, 2}}}),
        CBC(2, 2, {{0, 3}, {1, 3}, {2, 3}}, 12.75, 
                  {{{1, 3}, {2, 3}}, {{0, 3}, {2, 3}}, {{0, 3}, {1, 3}}}),
        CBC(2, 2, {{1, 2}, {1, 3}, {2, 3}}, 12.75, 
                  {{{1, 2}, {1, 3}}, {{1, 2}, {2, 3}}, {{1, 3}, {2, 3}}}),
        CBC(3, 1, {{0}, {0, 1}, {0, 2}, {0, 3}}, 6, 
                  {{{0}, {0, 2}, {0, 3}}, 
                   {{0}, {0, 1}, {0, 3}}, 
                   {{0}, {0, 1}, {0, 2}}, 
                   {{0, 1}, {0, 2}, {0, 3}}}),
        CBC(3, 1, {{0}, {0, 1}, {0, 2}, {1}, {1, 2}, {2}}, 12.5, 
                  {{{1}, {1, 2}, {2}}, 
                   {{0}, {0, 2}, {2}}, 
                   {{0}, {0, 1}, {1}}, 
                   {{0}, {0, 1}, {0, 2}}, 
                   {{0, 1}, {1}, {1, 2}}, 
                   {{0, 2}, {1, 2}, {2}}, 
                   {{0}, {1}, {2}}, 
                   {{0, 1}, {0, 2}, {1, 2}}}),
        CBC(3, 1, {{0}, {0, 1}, {0, 3}, {1}, {1, 3}, {3}}, 1.5, 
                  {{{1}, {1, 3}, {3}}, 
                   {{0}, {0, 3}, {3}}, 
                   {{0}, {0, 1}, {1}}, 
                   {{0}, {0, 1}, {0, 3}}, 
                   {{0, 1}, {1}, {1, 3}}, 
                   {{0, 3}, {1, 3}, {3}}, 
                   {{0}, {1}, {3}}, 
                   {{0, 1}, {0, 3}, {1, 3}}}),
        CBC(3, 1, {{0}, {0, 2}, {0, 3}, {2}, {2, 3}, {3}}, 1.5, 
                  {{{2}, {2, 3}, {3}}, 
                   {{0}, {0, 3}, {3}}, 
                   {{0}, {0, 2}, {2}}, 
                   {{0}, {0, 2}, {0, 3}}, 
                   {{0, 2}, {2}, {2, 3}}, 
                   {{0, 3}, {2, 3}, {3}}, 
                   {{0}, {2}, {3}}, 
                   {{0, 2}, {0, 3}, {2, 3}}}),
        CBC(3, 1, {{0}, {1}, {2}, {3}}, 12.75, 
                  {{{1}, {2}, {3}}, 
                   {{0}, {2}, {3}}, 
                   {{0}, {1}, {3}}, 
                   {{0}, {1}, {2}}}),
        CBC(3, 1, {{0, 1}, {1}, {1, 2}, {1, 3}}, 12.75, 
                  {{{1}, {1, 2}, {1, 3}}, 
                   {{0, 1}, {1}, {1, 3}}, 
                   {{0, 1}, {1}, {1, 2}}, 
                   {{0, 1}, {1, 2}, {1, 3}}}),
        CBC(3, 1, {{0, 2}, {1, 2}, {2}, {2, 3}}, 12.75, 
                  {{{1, 2}, {2}, {2, 3}}, 
                   {{0, 2}, {2}, {2, 3}}, 
                   {{0, 2}, {1, 2}, {2}}, 
                   {{0, 2}, {1, 2}, {2, 3}}}),
        CBC(3, 1, {{0, 3}, {1, 3}, {2, 3}, {3}}, 12.75, 
                  {{{1, 3}, {2, 3}, {3}}, 
                   {{0, 3}, {2, 3}, {3}}, 
                   {{0, 3}, {1, 3}, {3}}, 
                   {{0, 3}, {1, 3}, {2, 3}}}),
        CBC(3, 1, {{1}, {1, 2}, {1, 3}, {2}, {2, 3}, {3}}, 12.75, 
                  {{{2}, {2, 3}, {3}}, 
                   {{1}, {1, 3}, {3}}, 
                   {{1}, {1, 2}, {2}}, 
                   {{1}, {1, 2}, {1, 3}}, 
                   {{1, 2}, {2}, {2, 3}}, 
                   {{1, 3}, {2, 3}, {3}}, 
                   {{1}, {2}, {3}}, 
                   {{1, 2}, {1, 3}, {2, 3}}}),
        CBC(3, 2, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 12.75, 
                  {{{1, 2}, {1, 3}, {2, 3}}, 
                   {{0, 2}, {0, 3}, {2, 3}}, 
                   {{0, 1}, {0, 3}, {1, 3}}, 
                   {{0, 1}, {0, 2}, {1, 2}}, 
                   {{0, 1}, {0, 2}, {0, 3}}, 
                   {{0, 1}, {1, 2}, {1, 3}}, 
                   {{0, 2}, {1, 2}, {2, 3}}, 
                   {{0, 3}, {1, 3}, {2, 3}}}),
        CBC(4, 1, {{0}, {0, 1}, {0, 2}, {0, 3}, {1}, {1, 2}, {1, 3}, {2}, {2, 3}, {3}}, 12.75, 
                  {{{1}, {1, 2}, {1, 3}, {2}, {2, 3}, {3}}, 
                   {{0}, {0, 2}, {0, 3}, {2}, {2, 3}, {3}}, 
                   {{0}, {0, 1}, {0, 3}, {1}, {1, 3}, {3}}, 
                   {{0}, {0, 1}, {0, 2}, {1}, {1, 2}, {2}}, 
                   {{0}, {0, 1}, {0, 2}, {0, 3}}, 
                   {{0, 1}, {1}, {1, 2}, {1, 3}}, 
                   {{0, 2}, {1, 2}, {2}, {2, 3}}, 
                   {{0, 3}, {1, 3}, {2, 3}, {3}}, 
                   {{0}, {1}, {2}, {3}}, 
                   {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}}),
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 2);
    auto bf = rhomboidtiling.get_bifiltration();
    auto cm = rhomboidtiling.get_bifiltration_id_map();
    std::vector<CBC> cbifi;
    for (const auto& bc : bf) {
        cbifi.push_back(CBC(bc, cm));
    }
    std::sort(cbifi.begin(), cbifi.end());
    
    CHECK_THAT(cbifi_expected, UnorderedEquals(cbifi));
}


TEST_CASE("Half-integer slice critical rhomboid", "[halfint][minimal]") {

    Point p0(0, 0, 0);
    Point p1(0, 4, 4);
    Point p2(4, 4, 0);
    Point p3(4, 0, 4);

    std::vector<Point> points = {p0, p1, p2, p3};

    std::vector<std::vector<CVertex>> o15del_expected = {
        {{0, 1}, {0, 2}, {0, 3}, {1, 0}, {1, 2}, {1, 3}, 
         {2, 0}, {2, 1}, {2, 3}, {3, 0}, {3, 1}, {3, 2}}
    };

    std::vector<std::vector<CVertex>> o25del_expected = {
        {{0, 1, 2}, {0, 1, 3}, {0, 2, 1}, {0, 2, 3}, {0, 3, 1}, {0, 3, 2}, 
         {1, 2, 0}, {1, 2, 3}, {1, 3, 0}, {1, 3, 2}, {2, 3, 0}, {2, 3, 1}}
    };

    std::vector<std::vector<CVertex>> o35del_expected = {
        {{0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 3, 1}, {1, 2, 3, 0}}
    };

    auto rhomboidtiling = RhomboidTiling_3(points, 4);
    auto o15del = rhomboidtiling.get_halfint_slice_mosaic(1);
    auto o25del = rhomboidtiling.get_halfint_slice_mosaic(2);
    auto o35del = rhomboidtiling.get_halfint_slice_mosaic(3);
    CHECK_THAT(o15del, UnorderedEquals(o15del_expected));
    CHECK_THAT(o25del, UnorderedEquals(o25del_expected));
    CHECK_THAT(o35del, UnorderedEquals(o35del_expected));
}
