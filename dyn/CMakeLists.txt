add_library(dyn INTERFACE
        yoshida.h
        kahan.h
        halton.h
        newton.h
        circle.h
        verlet.h
        barnes_hut.h
)
target_include_directories(dyn INTERFACE .)
