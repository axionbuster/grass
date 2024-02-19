# Use WinMain() as entry point (if on Windows).
set(CMAKE_WIN32_EXECUTABLE TRUE)

add_executable(quadrantdemo main.cpp)

if(MSVC)
    target_compile_options(quadrantdemo PRIVATE /W4)
else()
    target_compile_options(quadrantdemo PRIVATE -Wall -Wextra -Wpedantic)
endif()

target_link_libraries(quadrantdemo dyn)
target_link_libraries(quadrantdemo raylib)

if (APPLE)
    target_link_libraries(grass "-framework IOKit")
    target_link_libraries(grass "-framework Cocoa")
    target_link_libraries(grass "-framework OpenGL")
endif ()
