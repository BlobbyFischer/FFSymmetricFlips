cmake_minimum_required(VERSION 4.2)
project(FFSymmetricFlips)

set(CMAKE_CXX_STANDARD 20)

add_executable(FFSymmetricFlips main.cpp
        F3Vec.h
        Symmetry336.h
        RankOneTensor.h
        Scheme.h
        FastRNG.h
        Utils.h)
