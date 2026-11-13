set(CMAKE_CXX_STANDARD 20)

add_library(JuniperSRC STATIC kernel.cpp kernel.h
                              SimData.cpp SimData.h
                              Tree.cpp Tree.h
                              ParticleOperators.cpp commons.h
                              GPUOperators.cpp GPUOperators.h
                              SimConfig.cpp SimConfig.h
                              ReportBuilders.cpp ReportBuilders.h
                              LeapfrogIntegrator.cpp LeapfrogIntegrator.h
                              toml.hpp)

find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found: ${OpenMP_FOUND}")
message(STATUS "OpenMP version: ${OpenMP_CXX_VERSION}")
message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}")

set(JUNIPER_GPU_ARCH "ccnative" CACHE STRING "nvc++ -gpu target (e.g. cc86)")
target_compile_options(JuniperSRC PUBLIC -mp=gpu -DTOML_HEADER_ONLY=0 -gpu=${JUNIPER_GPU_ARCH},nofma -Kieee -Mnofma -Minfo=mp -noswitcherror)
target_link_options(JuniperSRC PUBLIC -mp=gpu -DTOML_HEADER_ONLY=0 -gpu=${JUNIPER_GPU_ARCH} -noswitcherror)

target_include_directories(JuniperSRC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
