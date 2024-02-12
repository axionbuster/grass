include(FetchContent)
FetchContent_Declare(
        googletest
        # Specify the commit you depend on and update it regularly.
        URL https://github.com/google/googletest/archive/f8d7d77c06936315286eb55f8de22cd23c188571.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(GoogleTest)
include(GoogleTest)

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(units yoshida_test.cpp
        newton_test.cpp)
target_precompile_headers(units INTERFACE "gtest/gtest.h")
target_link_libraries(units gtest_main dyn)
gtest_discover_tests(units)
