cmake_minimum_required(VERSION 3.28.3)
project(Juniper2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LINKER_TYPE COMPILER)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(Juniper2 src/main.cpp)
