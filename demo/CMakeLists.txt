# Gravity

# Use WinMain() as entry point (if on Windows).
set(CMAKE_WIN32_EXECUTABLE TRUE)

add_executable(grass main.cpp
        Table.h
        user.h)

if (MSVC)
    target_compile_options(grass PRIVATE /W4)
    if (OPENMP)
        target_compile_options(grass PRIVATE /openmp:llvm)
    endif ()
else ()
    target_compile_options(grass PRIVATE -Wall -Wextra -Wpedantic)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND NOT APPLE)
        if (OPENMP)
            target_compile_options(grass PRIVATE -fopenmp=libiomp5)
            target_link_options(grass PRIVATE -fopenmp=libiomp5)
        endif ()
    else ()
        if (OPENMP)
            target_compile_options(grass PRIVATE -fopenmp)
            target_link_options(grass PRIVATE -fopenmp)
        endif ()
    endif ()
endif ()

target_link_libraries(grass raylib)
target_link_libraries(grass dyn)

# Header file too big; use precompiled header.
target_precompile_headers(grass INTERFACE <complex>)

if (APPLE)
    target_link_libraries(grass "-framework IOKit")
    target_link_libraries(grass "-framework Cocoa")
    target_link_libraries(grass "-framework OpenGL")
endif ()

if (EMSCRIPTEN)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --shell-file ${CMAKE_SOURCE_DIR}/demo/shell.html")
endif ()
