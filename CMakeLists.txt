cmake_minimum_required(VERSION 3.27)
project(grass)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(demo)
add_subdirectory(dyn)
add_subdirectory(quadrantdemo)

include(CTest)
if (BUILD_TESTING)
    add_subdirectory(tests)
endif ()
