cmake_minimum_required(VERSION 3.12)
project(Lab3)

set(CMAKE_CXX_STANDARD 14)

add_executable(Lab3 main.cpp HelperFunctions.cpp PRNG.h MultiplicativePRNG.cpp MultiplicativePRNG.h NormalDistributionModel.cpp NormalDistributionModel.h Cloneable.h LogNormalDistributionModel.cpp LogNormalDistributionModel.h CauchyDistributionModel.cpp CauchyDistributionModel.h)