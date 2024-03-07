# Gravity

# Use WinMain() as entry point (if on Windows).
set(CMAKE_WIN32_EXECUTABLE TRUE)

add_executable(hierarchydemo main.cpp "user.h")

if (MSVC)
    target_compile_options(hierarchydemo PRIVATE /W4)
else ()
    target_compile_options(hierarchydemo PRIVATE -Wall -Wextra -Wpedantic)
endif ()

target_link_libraries(hierarchydemo raylib)
target_link_libraries(hierarchydemo dyn)

if (APPLE)
    target_link_libraries(hierarchydemo "-framework IOKit")
    target_link_libraries(hierarchydemo "-framework Cocoa")
    target_link_libraries(hierarchydemo "-framework OpenGL")
endif ()
