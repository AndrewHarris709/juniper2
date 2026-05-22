enable_testing()

add_executable(tests
        test_blank.cpp
)

target_link_libraries(tests
        PRIVATE
        GTest::gtest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(tests)
