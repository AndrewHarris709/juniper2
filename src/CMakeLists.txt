set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LINKER_TYPE LLD)

add_library(JuniperSRC STATIC kernel.cpp kernel.h SimData.cpp SimData.h)

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found: ${OpenMP_FOUND}")
message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")

target_link_libraries(JuniperSRC PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(JuniperSRC PUBLIC -fopenmp -foffload=nvptx-none -fcf-protection=none -fno-stack-protector)
target_link_options(JuniperSRC PUBLIC -fopenmp -foffload=nvptx-none -fcf-protection=none -fno-stack-protector)

target_include_directories(JuniperSRC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
