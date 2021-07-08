cmake_minimum_required(VERSION 3.0.0)
project(PosTransTool CXX)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

file(GLOB SRCS ${CMAKE_SOURCE_DIR}/example/*.cpp)
file(GLOB HEAD ${CMAKE_SOURCE_DIR}/include/*.hpp)
# ==============================================================================
# Find Eigen
# ==============================================================================
find_package(Eigen3 REQUIRED)
if(NOT EIGEN3_FOUND)
	message(FATAL_ERROR "EIGEN3 not found.")
else()
	message("EIGEN3 library found.${EIGEN3_INCLUDE_DIR}")
endif()
set(EIGEN_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

add_executable(${PROJECT_NAME} ${SRCS} ${HEAD})
target_include_directories(${PROJECT_NAME} PUBLIC 
${EIGEN_INCLUDE_DIRS}
${CMAKE_SOURCE_DIR}/include)