set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LINKER_TYPE LLD)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0")
        message(FATAL_ERROR "This project requires GCC 16 or newer. (found GCC ${CMAKE_CXX_COMPILER_VERSION})")
    endif()
endif()

add_library(JuniperSRC STATIC kernel.cpp kernel.h
                              SimData.cpp SimData.h
                              Tree.cpp Tree.h
                              ParticleOperators.cpp commons.h
                              GPUOperators.cpp GPUOperators.h
                              SimConfig.cpp SimConfig.h
                              ReportBuilders.cpp ReportBuilders.h
                              toml.hpp)

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found: ${OpenMP_FOUND}")
message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")

target_link_libraries(JuniperSRC PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(JuniperSRC PUBLIC -fopenmp -foffload=nvptx-none -fno-fast-math -ffp-contract=off -fno-associative-math -fcf-protection=none -fno-stack-protector -foffload-options=nvptx-none="-lm")
target_link_options(JuniperSRC PUBLIC -fopenmp -foffload=nvptx-none -fno-fast-math -ffp-contract=off -fno-associative-math -fcf-protection=none -fno-stack-protector -foffload-options=nvptx-none="-lm")

target_include_directories(JuniperSRC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
