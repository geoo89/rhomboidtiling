/*
 * Copyright (c) 2019 Georg Osang
 * Distributed under the MIT License, see LICENCE.md
 */

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "dimensional_traits_2.h"
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

#include "testutils.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel                    K;
typedef K::FT                                                               FT;
typedef DimensionalTraits_2<K>                                             Dt2;
typedef RhomboidTiling<Dt2>                                   RhomboidTiling_2;
typedef RhomboidTiling_2::Point                                          Point;
typedef CombinatorialBifiltrationCell<FT>                                  CBC;


TEST_CASE("Standard toy example for orderk", "[orderk][toy]") {

    std::vector<Point> points = {
        Point(156.006,705.854),
        Point(215.257,732.63),
        Point(283.108,707.272),
        Point(244.042,670.948),
        Point(366.035,687.396),
        Point(331.768,625.715),
        Point(337.936,559.92),
        Point(249.525,582.537),
        Point(187.638,556.13),
        Point(165.912,631.197)
    };

    std::vector<std::vector<CVertex>> o2del_expected = {
        {{0, 1}, {0, 3}, {0, 9}},
        {{0, 1}, {0, 3}, {1, 3}},
        {{0, 1}, {1, 2}, {1, 3}},
        {{0, 3}, {0, 9}, {3, 9}},
        {{0, 3}, {1, 3}, {3, 9}},
        {{0, 9}, {3, 9}, {8, 9}},
        {{1, 2}, {1, 3}, {2, 3}},
        {{1, 2}, {1, 4}, {2, 4}},
        {{1, 2}, {2, 3}, {2, 4}},
        {{1, 3}, {2, 3}, {3, 9}},
        {{2, 3}, {2, 4}, {2, 5}},
        {{2, 3}, {2, 5}, {3, 5}},
        {{2, 3}, {3, 5}, {3, 7}},
        {{2, 3}, {3, 7}, {3, 9}},
        {{2, 4}, {2, 5}, {4, 5}},
        {{2, 5}, {3, 5}, {4, 5}},
        {{3, 5}, {3, 7}, {5, 7}},
        {{3, 5}, {4, 5}, {5, 6}},
        {{3, 5}, {5, 6}, {5, 7}},
        {{3, 7}, {3, 9}, {7, 9}},
        {{3, 7}, {5, 7}, {7, 8}},
        {{3, 7}, {7, 8}, {7, 9}},
        {{3, 9}, {7, 9}, {8, 9}},
        {{4, 5}, {4, 6}, {5, 6}},
        {{5, 6}, {5, 7}, {6, 7}},
        {{5, 7}, {6, 7}, {7, 8}},
        {{6, 7}, {6, 8}, {7, 8}},
        {{7, 8}, {7, 9}, {8, 9}}
    };

    std::vector<std::vector<CVertex>> o4del_expected = {
        {{0, 1, 2, 3}, {0, 1, 2, 4}, {1, 2, 3, 4}},
        {{0, 1, 2, 3}, {0, 1, 2, 9}, {0, 1, 3, 9}},
        {{0, 1, 2, 3}, {0, 1, 3, 9}, {1, 2, 3, 9}},
        {{0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 5}},
        {{0, 1, 2, 3}, {1, 2, 3, 5}, {1, 2, 3, 7}},
        {{0, 1, 2, 3}, {1, 2, 3, 7}, {1, 2, 3, 9}},
        {{0, 1, 3, 9}, {0, 1, 8, 9}, {0, 3, 8, 9}},
        {{0, 1, 3, 9}, {0, 3, 7, 9}, {0, 3, 8, 9}},
        {{0, 1, 3, 9}, {0, 3, 7, 9}, {1, 3, 7, 9}},
        {{0, 1, 3, 9}, {1, 2, 3, 9}, {1, 3, 7, 9}},
        {{0, 1, 8, 9}, {0, 3, 8, 9}, {0, 7, 8, 9}},
        {{0, 3, 7, 9}, {0, 3, 8, 9}, {3, 7, 8, 9}},
        {{0, 3, 7, 9}, {1, 3, 7, 9}, {2, 3, 7, 9}},
        {{0, 3, 7, 9}, {2, 3, 7, 9}, {3, 7, 8, 9}},
        {{0, 3, 8, 9}, {0, 7, 8, 9}, {3, 7, 8, 9}},
        {{0, 7, 8, 9}, {3, 7, 8, 9}, {6, 7, 8, 9}},
        {{1, 2, 3, 4}, {1, 2, 3, 5}, {2, 3, 4, 5}},
        {{1, 2, 3, 4}, {1, 2, 4, 5}, {2, 3, 4, 5}},
        {{1, 2, 3, 5}, {1, 2, 3, 7}, {2, 3, 5, 7}},
        {{1, 2, 3, 5}, {2, 3, 4, 5}, {2, 3, 5, 7}},
        {{1, 2, 3, 7}, {1, 2, 3, 9}, {2, 3, 7, 9}},
        {{1, 2, 3, 7}, {2, 3, 5, 7}, {2, 3, 7, 9}},
        {{1, 2, 3, 9}, {1, 3, 7, 9}, {2, 3, 7, 9}},
        {{1, 2, 4, 5}, {2, 3, 4, 5}, {2, 4, 5, 6}},
        {{2, 3, 4, 5}, {2, 3, 5, 7}, {3, 4, 5, 7}},
        {{2, 3, 4, 5}, {2, 4, 5, 6}, {3, 4, 5, 6}},
        {{2, 3, 4, 5}, {3, 4, 5, 6}, {3, 4, 5, 7}},
        {{2, 3, 5, 7}, {2, 3, 7, 9}, {3, 5, 7, 9}},
        {{2, 3, 5, 7}, {3, 4, 5, 7}, {3, 5, 6, 7}},
        {{2, 3, 5, 7}, {3, 5, 6, 7}, {3, 5, 7, 8}},
        {{2, 3, 5, 7}, {3, 5, 7, 8}, {3, 5, 7, 9}},
        {{2, 3, 7, 9}, {3, 5, 7, 9}, {3, 7, 8, 9}},
        {{2, 4, 5, 6}, {3, 4, 5, 6}, {4, 5, 6, 7}},
        {{3, 4, 5, 6}, {3, 4, 5, 7}, {3, 5, 6, 7}},
        {{3, 4, 5, 6}, {3, 5, 6, 7}, {4, 5, 6, 7}},
        {{3, 5, 6, 7}, {3, 5, 7, 8}, {5, 6, 7, 8}},
        {{3, 5, 6, 7}, {4, 5, 6, 7}, {5, 6, 7, 8}},
        {{3, 5, 7, 8}, {3, 5, 7, 9}, {3, 7, 8, 9}},
        {{3, 5, 7, 8}, {3, 6, 7, 8}, {3, 7, 8, 9}},
        {{3, 5, 7, 8}, {3, 6, 7, 8}, {5, 6, 7, 8}},
        {{3, 6, 7, 8}, {3, 7, 8, 9}, {6, 7, 8, 9}},
        {{3, 6, 7, 8}, {5, 6, 7, 8}, {6, 7, 8, 9}},
    };

    std::vector<std::vector<CVertex>> o7del_expected = {
        {{0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 4, 5, 7}, {0, 1, 2, 3, 4, 5, 9}},
        {{0, 1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 4, 5, 7}, {1, 2, 3, 4, 5, 6, 7}},
        {{0, 1, 2, 3, 4, 5, 7}, {0, 1, 2, 3, 4, 5, 9}, {0, 1, 2, 3, 5, 7, 9}},
        {{0, 1, 2, 3, 4, 5, 7}, {0, 1, 2, 3, 5, 7, 9}, {1, 2, 3, 4, 5, 7, 9}},
        {{0, 1, 2, 3, 4, 5, 7}, {1, 2, 3, 4, 5, 6, 7}, {1, 2, 3, 4, 5, 7, 9}},
        {{0, 1, 2, 3, 4, 5, 9}, {0, 1, 2, 3, 4, 7, 9}, {0, 1, 2, 3, 4, 8, 9}},
        {{0, 1, 2, 3, 4, 5, 9}, {0, 1, 2, 3, 4, 7, 9}, {0, 1, 2, 3, 5, 7, 9}},
        {{0, 1, 2, 3, 4, 7, 9}, {0, 1, 2, 3, 4, 8, 9}, {0, 1, 2, 3, 7, 8, 9}},
        {{0, 1, 2, 3, 4, 7, 9}, {0, 1, 2, 3, 5, 7, 9}, {0, 1, 2, 3, 7, 8, 9}},
        {{0, 1, 2, 3, 5, 7, 9}, {0, 1, 2, 3, 7, 8, 9}, {1, 2, 3, 5, 7, 8, 9}},
        {{0, 1, 2, 3, 5, 7, 9}, {1, 2, 3, 4, 5, 7, 9}, {1, 2, 3, 5, 7, 8, 9}},
        {{0, 1, 2, 3, 7, 8, 9}, {0, 1, 3, 5, 7, 8, 9}, {0, 1, 3, 6, 7, 8, 9}},
        {{0, 1, 2, 3, 7, 8, 9}, {0, 1, 3, 5, 7, 8, 9}, {0, 2, 3, 5, 7, 8, 9}},
        {{0, 1, 2, 3, 7, 8, 9}, {0, 2, 3, 5, 7, 8, 9}, {1, 2, 3, 5, 7, 8, 9}},
        {{0, 1, 3, 5, 7, 8, 9}, {0, 1, 3, 6, 7, 8, 9}, {0, 3, 5, 6, 7, 8, 9}},
        {{0, 1, 3, 5, 7, 8, 9}, {0, 2, 3, 5, 7, 8, 9}, {0, 3, 5, 6, 7, 8, 9}},
        {{0, 2, 3, 5, 7, 8, 9}, {0, 3, 5, 6, 7, 8, 9}, {2, 3, 5, 6, 7, 8, 9}},
        {{0, 2, 3, 5, 7, 8, 9}, {1, 2, 3, 5, 7, 8, 9}, {2, 3, 5, 6, 7, 8, 9}},
        {{0, 3, 5, 6, 7, 8, 9}, {2, 3, 5, 6, 7, 8, 9}, {3, 4, 5, 6, 7, 8, 9}},
        {{1, 2, 3, 4, 5, 6, 7}, {1, 2, 3, 4, 5, 7, 9}, {2, 3, 4, 5, 6, 7, 9}},
        {{1, 2, 3, 4, 5, 6, 7}, {2, 3, 4, 5, 6, 7, 8}, {2, 3, 4, 5, 6, 7, 9}},
        {{1, 2, 3, 4, 5, 7, 9}, {1, 2, 3, 5, 6, 7, 9}, {1, 2, 3, 5, 7, 8, 9}},
        {{1, 2, 3, 4, 5, 7, 9}, {1, 2, 3, 5, 6, 7, 9}, {2, 3, 4, 5, 6, 7, 9}},
        {{1, 2, 3, 5, 6, 7, 9}, {1, 2, 3, 5, 7, 8, 9}, {2, 3, 5, 6, 7, 8, 9}},
        {{1, 2, 3, 5, 6, 7, 9}, {2, 3, 4, 5, 6, 7, 9}, {2, 3, 5, 6, 7, 8, 9}},
        {{2, 3, 4, 5, 6, 7, 8}, {2, 3, 4, 5, 6, 7, 9}, {2, 3, 5, 6, 7, 8, 9}},
        {{2, 3, 4, 5, 6, 7, 8}, {2, 3, 5, 6, 7, 8, 9}, {3, 4, 5, 6, 7, 8, 9}},
    };

    std::vector<std::vector<CVertex>> o9del_expected = {
        {{0, 1, 2, 3, 4, 5, 6, 7, 8}, {0, 1, 2, 3, 4, 5, 6, 7, 9}, {1, 2, 3, 4, 5, 6, 7, 8, 9}},
        {{0, 1, 2, 3, 4, 5, 6, 7, 9}, {0, 1, 2, 3, 4, 5, 7, 8, 9}, {1, 2, 3, 4, 5, 6, 7, 8, 9}},
        {{0, 1, 2, 3, 4, 5, 7, 8, 9}, {0, 1, 2, 3, 5, 6, 7, 8, 9}, {1, 2, 3, 4, 5, 6, 7, 8, 9}},
        {{0, 1, 2, 3, 5, 6, 7, 8, 9}, {0, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 2, 3, 4, 5, 6, 7, 8, 9}}
    };

    auto rhomboidtiling = RhomboidTiling_2(points, 10);
    auto o2del = rhomboidtiling.get_slice_mosaic(2);
    auto o4del = rhomboidtiling.get_slice_mosaic(4);
    auto o7del = rhomboidtiling.get_slice_mosaic(7);
    auto o9del = rhomboidtiling.get_slice_mosaic(9);
    CHECK_THAT(o2del, UnorderedEquals(o2del_expected));
    CHECK_THAT(o4del, UnorderedEquals(o4del_expected));
    CHECK_THAT(o7del, UnorderedEquals(o7del_expected));
    CHECK_THAT(o9del, UnorderedEquals(o9del_expected));
}


TEST_CASE("Radius for small example", "[orderkfilt]") {

    std::vector<Point> points = {
        Point(0,1),
        Point(0,3),
        Point(1,2),
        Point(3,3),
        Point(4,0),
        Point(4,2),
    };

    std::vector<CellWithRadius<FT>> o1del_expected = {
        CellWithRadius<FT>({{1}}, 0),
        CellWithRadius<FT>({{2}}, 0),
        CellWithRadius<FT>({{3}}, 0),
        CellWithRadius<FT>({{5}}, 0),
        CellWithRadius<FT>({{4}}, 0),
        CellWithRadius<FT>({{0}}, 0),
        CellWithRadius<FT>({{3}, {5}}, 0.5),
        CellWithRadius<FT>({{0}, {2}}, 0.5),
        CellWithRadius<FT>({{1}, {2}}, 0.5),
        CellWithRadius<FT>({{0}, {1}}, 1),
        CellWithRadius<FT>({{0}, {1}, {2}}, 1),
        CellWithRadius<FT>({{4}, {5}}, 1),
        CellWithRadius<FT>({{2}, {3}}, 1.25),
        CellWithRadius<FT>({{2}, {5}}, 2.5),
        CellWithRadius<FT>({{2}, {3}, {5}}, 2.5),
        CellWithRadius<FT>({{1}, {3}}, 2.5),
        CellWithRadius<FT>({{1}, {2}, {3}}, 2.5),
        CellWithRadius<FT>({{2}, {4}}, 3.25),
        CellWithRadius<FT>({{2}, {4}, {5}}, 3.25),
        CellWithRadius<FT>({{0}, {4}}, 4.42),
        CellWithRadius<FT>({{0}, {2}, {4}}, 4.42),
    };

    std::vector<CellWithRadius<FT>> o2del_expected = {
        CellWithRadius<FT>({{0, 2}}, 0.5),
        CellWithRadius<FT>({{1, 2}}, 0.5),
        CellWithRadius<FT>({{3, 5}}, 0.5),
        CellWithRadius<FT>({{4, 5}}, 1),
        CellWithRadius<FT>({{0, 1}}, 1),
        CellWithRadius<FT>({{0, 1}, {1, 2}}, 1),
        CellWithRadius<FT>({{0, 1}, {0, 2}}, 1),
        CellWithRadius<FT>({{0, 2}, {1, 2}}, 1),
        CellWithRadius<FT>({{0, 1}, {0, 2}, {1, 2}}, 1),
        CellWithRadius<FT>({{2, 3}}, 1.25),
        CellWithRadius<FT>({{1, 2}, {2, 3}}, 2.25),
        CellWithRadius<FT>({{2, 3}, {3, 5}}, 2.25),
        CellWithRadius<FT>({{3, 5}, {4, 5}}, 2.5),
        CellWithRadius<FT>({{2, 5}}, 2.5),
        CellWithRadius<FT>({{2, 3}, {2, 5}}, 2.5),
        CellWithRadius<FT>({{2, 5}, {3, 5}}, 2.5),
        CellWithRadius<FT>({{2, 3}, {2, 5}, {3, 5}}, 2.5),
        CellWithRadius<FT>({{1, 3}}, 2.5),
        CellWithRadius<FT>({{1, 2}, {1, 3}}, 2.5),
        CellWithRadius<FT>({{1, 3}, {2, 3}}, 2.5),
        CellWithRadius<FT>({{1, 2}, {1, 3}, {2, 3}}, 2.5),
        CellWithRadius<FT>({{2, 4}}, 3.25),
        CellWithRadius<FT>({{2, 4}, {2, 5}}, 3.25),
        CellWithRadius<FT>({{2, 4}, {4, 5}}, 3.25),
        CellWithRadius<FT>({{2, 5}, {4, 5}}, 3.25),
        CellWithRadius<FT>({{2, 4}, {2, 5}, {4, 5}}, 3.25),
        CellWithRadius<FT>({{0, 2}, {2, 3}}, 3.25),
        CellWithRadius<FT>({{0, 2}, {1, 2}, {2, 3}}, 3.25),
        CellWithRadius<FT>({{2, 5}, {3, 5}, {4, 5}}, FT(162.5)/49),
        CellWithRadius<FT>({{0, 2}, {2, 4}}, 4.25),
        CellWithRadius<FT>({{0, 4}}, 4.42),
        CellWithRadius<FT>({{0, 2}, {0, 4}}, 4.42),
        CellWithRadius<FT>({{0, 4}, {2, 4}}, 4.42),
        CellWithRadius<FT>({{0, 2}, {0, 4}, {2, 4}}, 4.42),
        CellWithRadius<FT>({{0, 2}, {2, 5}}, 4.42),
        CellWithRadius<FT>({{0, 2}, {2, 3}, {2, 5}}, 4.42),
        CellWithRadius<FT>({{0, 2}, {2, 4}, {2, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 4}, {4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 4}, {2, 4}, {4, 5}}, 8.5),
        CellWithRadius<FT>({{1, 3}, {3, 5}}, 8.5),
        CellWithRadius<FT>({{1, 3}, {2, 3}, {3, 5}}, 8.5),
    };

    std::vector<CellWithRadius<FT>> o3del_expected = {
        CellWithRadius<FT>({{0, 1, 2}}, 1),
        CellWithRadius<FT>({{2, 3, 5}}, 2.25),
        CellWithRadius<FT>({{1, 2, 3}}, 2.25),
        CellWithRadius<FT>({{3, 4, 5}}, 2.5),
        CellWithRadius<FT>({{0, 2, 3}}, 3.25),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}}, 3.25),
        CellWithRadius<FT>({{0, 2, 3}, {1, 2, 3}}, 3.25),
        CellWithRadius<FT>({{0, 1, 2}, {1, 2, 3}}, 3.25),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}, {1, 2, 3}}, 3.25),
        CellWithRadius<FT>({{2, 4, 5}}, 3.25),
        CellWithRadius<FT>({{2, 3, 5}, {3, 4, 5}}, FT(162.5)/49),
        CellWithRadius<FT>({{2, 4, 5}, {3, 4, 5}}, FT(162.5)/49),
        CellWithRadius<FT>({{2, 3, 5}, {2, 4, 5}}, FT(162.5)/49),
        CellWithRadius<FT>({{2, 3, 5}, {2, 4, 5}, {3, 4, 5}}, FT(162.5)/49),
        CellWithRadius<FT>({{0, 2, 3}, {2, 3, 5}}, 4.25),
        CellWithRadius<FT>({{1, 2, 3}, {2, 3, 5}}, 4.25),
        CellWithRadius<FT>({{0, 2, 4}}, 4.25),
        CellWithRadius<FT>({{0, 2, 5}}, 4.42),
        CellWithRadius<FT>({{0, 2, 5}, {2, 3, 5}}, 4.42),
        CellWithRadius<FT>({{0, 2, 3}, {0, 2, 5}}, 4.42),
        CellWithRadius<FT>({{0, 2, 3}, {0, 2, 5}, {2, 3, 5}}, 4.42),
        CellWithRadius<FT>({{0, 2, 3}, {1, 2, 3}, {2, 3, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 4}, {2, 4, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 5}, {2, 4, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 4}, {0, 2, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 4}, {0, 2, 5}, {2, 4, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 5}, {2, 3, 5}, {2, 4, 5}}, FT(552.5)/121),
        CellWithRadius<FT>({{0, 2, 3}, {0, 2, 4}}, 5),
        CellWithRadius<FT>({{0, 2, 3}, {0, 2, 4}, {0, 2, 5}}, 5),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 4}}, FT(62.5)/9),
        CellWithRadius<FT>({{0, 1, 2}, {0, 2, 3}, {0, 2, 4}}, FT(62.5)/9),
        CellWithRadius<FT>({{0, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 4, 5}, {2, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 2, 4}, {0, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 2, 4}, {0, 4, 5}, {2, 4, 5}}, 8.5),
        CellWithRadius<FT>({{1, 3, 5}}, 8.5),
        CellWithRadius<FT>({{1, 2, 3}, {1, 3, 5}}, 8.5),
        CellWithRadius<FT>({{1, 3, 5}, {2, 3, 5}}, 8.5),
        CellWithRadius<FT>({{1, 2, 3}, {1, 3, 5}, {2, 3, 5}}, 8.5),
        CellWithRadius<FT>({{0, 4, 5}, {3, 4, 5}}, 32.5),
        CellWithRadius<FT>({{0, 4, 5}, {2, 4, 5}, {3, 4, 5}}, 32.5),
        CellWithRadius<FT>({{1, 3, 5}, {3, 4, 5}}, 162.5),
        CellWithRadius<FT>({{1, 3, 5}, {2, 3, 5}, {3, 4, 5}}, 162.5),
    };

    std::vector<CellWithRadius<FT>> o4del_expected = {
        CellWithRadius<FT>({{0, 1, 2, 3}}, 3.25),
        CellWithRadius<FT>({{2, 3, 4, 5}}, FT(162.5)/49),
        CellWithRadius<FT>({{1, 2, 3, 5}}, 4.25),
        CellWithRadius<FT>({{0, 2, 3, 5}}, 4.25),
        CellWithRadius<FT>({{0, 1, 2, 3}, {1, 2, 3, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 3, 5}, {1, 2, 3, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 2, 3, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 2, 3, 5}, {1, 2, 3, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 4, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 4, 5}, {2, 3, 4, 5}}, FT(552.5)/121),
        CellWithRadius<FT>({{0, 2, 3, 5}, {0, 2, 4, 5}}, FT(552.5)/121),
        CellWithRadius<FT>({{0, 2, 3, 5}, {2, 3, 4, 5}}, FT(552.5)/121),
        CellWithRadius<FT>({{0, 2, 3, 5}, {0, 2, 4, 5}, {2, 3, 4, 5}}, FT(552.5)/121),
        CellWithRadius<FT>({{0, 2, 3, 4}}, 5),
        CellWithRadius<FT>({{0, 2, 3, 4}, {0, 2, 4, 5}}, 5),
        CellWithRadius<FT>({{0, 2, 3, 4}, {0, 2, 3, 5}}, 5),
        CellWithRadius<FT>({{0, 2, 3, 4}, {0, 2, 3, 5}, {0, 2, 4, 5}}, 5),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 2, 3, 4}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 2, 3, 4}, {0, 2, 3, 5}}, 6.640625),
        CellWithRadius<FT>({{1, 2, 3, 5}, {2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 2, 3, 5}, {1, 2, 3, 5}, {2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 4}}, FT(62.5)/9),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}}, FT(62.5)/9),
        CellWithRadius<FT>({{0, 1, 2, 4}, {0, 2, 3, 4}}, FT(62.5)/9),
        CellWithRadius<FT>({{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 2, 3, 4}}, FT(62.5)/9),
        CellWithRadius<FT>({{0, 1, 2, 4}, {0, 2, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 1, 2, 4}, {0, 2, 3, 4}, {0, 2, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 3, 4, 5}}, 32.5),
        CellWithRadius<FT>({{0, 2, 4, 5}, {0, 3, 4, 5}}, 32.5),
        CellWithRadius<FT>({{0, 3, 4, 5}, {2, 3, 4, 5}}, 32.5),
        CellWithRadius<FT>({{0, 2, 4, 5}, {0, 3, 4, 5}, {2, 3, 4, 5}}, 32.5),
        CellWithRadius<FT>({{1, 3, 4, 5}}, 162.5),
        CellWithRadius<FT>({{1, 2, 3, 5}, {1, 3, 4, 5}}, 162.5),
        CellWithRadius<FT>({{1, 3, 4, 5}, {2, 3, 4, 5}}, 162.5),
        CellWithRadius<FT>({{1, 2, 3, 5}, {1, 3, 4, 5}, {2, 3, 4, 5}}, 162.5),
    };

    std::vector<CellWithRadius<FT>> o5del_expected = {
        CellWithRadius<FT>({{0, 1, 2, 3, 5}}, 4.515625),
        CellWithRadius<FT>({{0, 2, 3, 4, 5}}, FT(552.5)/121),
        CellWithRadius<FT>({{0, 1, 2, 3, 5}, {0, 2, 3, 4, 5}}, 6.25),
        CellWithRadius<FT>({{1, 2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 2, 3, 4, 5}, {1, 2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 3, 5}, {1, 2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 3, 5}, {0, 2, 3, 4, 5}, {1, 2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 3, 4}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 3, 4}, {0, 2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 3, 4}, {0, 1, 2, 3, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 3, 4}, {0, 1, 2, 3, 5}, {0, 2, 3, 4, 5}}, 6.640625),
        CellWithRadius<FT>({{0, 1, 2, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 1, 2, 4, 5}, {0, 2, 3, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 1, 2, 3, 4}, {0, 1, 2, 4, 5}}, 8.5),
        CellWithRadius<FT>({{0, 1, 2, 3, 4}, {0, 1, 2, 4, 5}, {0, 2, 3, 4, 5}}, 8.5),
    };

    auto rhomboidtiling = RhomboidTiling_2(points, 6);
    auto o1del = rhomboidtiling.get_slice_filtration(1);
    auto o2del = rhomboidtiling.get_slice_filtration(2);
    auto o3del = rhomboidtiling.get_slice_filtration(3);
    auto o4del = rhomboidtiling.get_slice_filtration(4);
    auto o5del = rhomboidtiling.get_slice_filtration(5);
    CHECK_THAT(o1del, UnorderedEquals(o1del_expected));
    CHECK_THAT(o2del, UnorderedEquals(o2del_expected));
    CHECK_THAT(o3del, UnorderedEquals(o3del_expected));
    CHECK_THAT(o4del, UnorderedEquals(o4del_expected));
    CHECK_THAT(o5del, UnorderedEquals(o5del_expected));
}
