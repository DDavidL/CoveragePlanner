cmake_minimum_required(VERSION 3.5)
project(CoveragePlanner)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3 3.4.0 REQUIRED)
message(STATUS "Please ensure Eigen version 3.4.0 headers are installed and accessible in your build environment.")
find_package(OpenCV 4.10.0 REQUIRED)
message(STATUS "Please ensure OpenCV version 4.10.0 is installed and accessible in your build environment.")
find_package(CGAL 6.0.1 REQUIRED COMPONENTS Core)
message(STATUS "Please ensure CGAL version 6.0.1 and its dependencies (GMP, MPFR, Eigen 3.4.0) are installed and accessible in your build environment.")

include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})

file(GLOB_RECURSE srcs "src/*.cc" "src/*.cpp")
file(GLOB_RECURSE hdrs "include/*.h")

message("find source files: ${srcs}")
message("find headers: ${hdrs}")

include_directories(include)

add_executable(CoveragePlanner ${srcs} ${hdrs})

target_link_libraries(CoveragePlanner
        ${OpenCV_LIBS}
        CGAL::CGAL
        CGAL::CGAL_Core
        )

