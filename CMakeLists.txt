cmake_minimum_required(VERSION 3.0.0)
project(pos2ecf CXX)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11")
set(CMAKE_VERBOSE_MAKEFILE ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

aux_source_directory(./src SRC_DIR)

#Eigen3
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
  set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
endif()

# set MD to MT
if(OPTION_BUILD_MT)
    set(CompilerFlags
            CMAKE_CXX_FLAGS
            CMAKE_CXX_FLAGS_DEBUG
            CMAKE_CXX_FLAGS_RELEASE
            CMAKE_C_FLAGS
            CMAKE_C_FLAGS_DEBUG
            CMAKE_C_FLAGS_RELEASE
            )
foreach(CompilerFlags ${CompilerFlags})
    string(REPLACE "/MD" "/MT" ${CompilerFlags} "${${CompilerFlags}}")
    string(REPLACE "/Od" "/O2" ${CompilerFlags} "${${CompilerFlags}}")
endforeach()
endif()

include_directories(${EIGEN3_INCLUDE_DIR})
add_executable(${PROJECT_NAME} ${SRC_DIR})