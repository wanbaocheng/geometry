cmake_minimum_required(VERSION 3.12)
project(geometry LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(demo examples/demo.cpp)
target_link_libraries(demo)
