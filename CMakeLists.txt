cmake_minimum_required(VERSION 3.10)
project(CoveragePlanner)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenCV 4 REQUIRED COMPONENTS core imgproc highgui)
find_package(CGAL REQUIRED COMPONENTS Core)



file(GLOB_RECURSE srcs "src/*.cc" "src/*.cpp")
file(GLOB_RECURSE hdrs "include/*.h")


message("find source files: ${srcs}")
message("find headers: ${hdrs}")

add_executable(CoveragePlanner ${srcs} ${hdrs})

target_include_directories(CoveragePlanner PRIVATE include)
target_link_libraries(CoveragePlanner PRIVATE
        Eigen3::Eigen
        ${OpenCV_LIBS}
        CGAL::CGAL
        CGAL::CGAL_Core
        )

enable_testing()
add_subdirectory(tests)

