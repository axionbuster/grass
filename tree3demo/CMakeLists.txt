# Use WinMain() as entry point (if on Windows).
set(CMAKE_WIN32_EXECUTABLE TRUE)

add_executable(tree3demo main.cpp)

if (MSVC)
    target_compile_options(tree3demo PRIVATE /W4)
else ()
    target_compile_options(tree3demo PRIVATE -Wall -Wextra -Wpedantic)
endif ()

target_link_libraries(tree3demo dyn)
target_link_libraries(tree3demo raylib)

if (APPLE)
    target_link_libraries(grass "-framework IOKit")
    target_link_libraries(grass "-framework Cocoa")
    target_link_libraries(grass "-framework OpenGL")
endif ()
