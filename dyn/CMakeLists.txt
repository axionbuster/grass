cmake_minimum_required(VERSION 3.27)

add_library(dyn INTERFACE
        yoshida.h
        kahan.h
        halton.h
        newton.h
        circle.h
        verlet.h)
target_include_directories(dyn INTERFACE .)
