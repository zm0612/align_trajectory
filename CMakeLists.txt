cmake_minimum_required(VERSION 2.8)
project(align_trajectory)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE "release")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)

find_package(Eigen3)
find_package(Pangolin REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${Pangolin_INCLUDE_DIRS})
include_directories("./include")

add_executable(align_trajectory main.cpp)
target_link_libraries(align_trajectory
        ${Pangolin_LIBRARIES}
)
