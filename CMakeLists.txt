cmake_minimum_required(VERSION 3.28.3)
project(Juniper2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LINKER_TYPE COMPILER)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(Juniper2 src/main.cpp)
target_link_libraries(Juniper2 PUBLIC JuniperSRC)
