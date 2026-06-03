include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(IOTest test_io.cpp)
target_link_libraries(IOTest
        PRIVATE
        GTest::gtest_main
        GTest::gmock_main
        JuniperSRC
)
add_test(NAME IO_Tests COMMAND IOTest)

include(GoogleTest)
gtest_discover_tests(IOTest)
