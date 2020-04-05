# CMake for order-k Delaunay triangulations using CGAL.
# Includes a commandline tool (main) and tests (test).

project (rhomboid)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_definitions(-std=c++11)

find_package(CGAL COMPONENTS)
include(${CGAL_USE_FILE}) 

if ( CGAL_FOUND )

set(DEPENDENCIES src/rhomboid.cpp src/utils.cpp)

add_executable(main src/main.cpp ${DEPENDENCIES})
add_executable(tests_2d src/tests_2d.cpp ${DEPENDENCIES})
add_executable(tests_3d src/tests_3d.cpp ${DEPENDENCIES})
add_to_cached_list( CGAL_EXECUTABLE_TARGETS main)
target_link_libraries( main ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL, and will not be compiled.")

endif()