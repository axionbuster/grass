# CMakeList.txt : CMake project for grass, include source and define
# project specific logic here.
#

add_library(wow)

target_sources(wow
    PUBLIC
        FILE_SET CXX_MODULES FILES
            wow.cppm
    PRIVATE
        wow.cpp
)

# Add source to this project's executable.
add_executable (grass grass.cpp)

# Link.
target_link_libraries(grass wow)

# TODO: Add tests and install targets if needed.
