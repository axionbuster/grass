# Gravity

# Use WinMain() as entry point (if on Windows).
set(CMAKE_WIN32_EXECUTABLE TRUE)

add_executable(grass main.cpp
        Table.h
        user.h)

if (MSVC)
    target_compile_options(grass PRIVATE /W4)
else ()
    target_compile_options(grass PRIVATE -Wall -Wextra -Wpedantic)
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
