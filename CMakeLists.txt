cmake_minimum_required (VERSION 3.1)
project (lambdatwist)

if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()


set (CMAKE_CXX_STANDARD 14)

add_compile_options(-Wall -Werror -Wno-sign-compare -Wno-unused-variable)

find_package (Eigen3 REQUIRED NO_MODULE)

add_library(lambdatwist STATIC lambdatwist/p3p.cc)
target_link_libraries (lambdatwist Eigen3::Eigen)

add_executable (lambdatwist_test lambdatwist/p3p_test.cc)
target_link_libraries (lambdatwist_test LINK_PUBLIC lambdatwist Eigen3::Eigen)
