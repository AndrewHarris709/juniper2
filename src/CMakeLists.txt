set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LINKER_TYPE LLD)

add_executable(my_app main.cpp kernel.cpp)
target_include_directories(my_app PUBLIC .)

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found: ${OpenMP_FOUND}")
message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")

target_link_libraries(my_app PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(my_app PUBLIC -fopenmp -foffload=nvptx-none -fcf-protection=none -fno-stack-protector)
target_link_options(my_app PUBLIC -fopenmp -foffload=nvptx-none -fcf-protection=none -fno-stack-protector)
