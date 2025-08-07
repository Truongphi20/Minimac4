cmake_minimum_required(VERSION 3.14)
project(minimac4 VERSION 4.1.6)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type (Debug|Release|RelWithDebInfo|MinSizeRel)" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)

add_definitions(-DVERSION="${PROJECT_VERSION}")

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

## Install OMP
include(FetchContent)

FetchContent_Declare(
  omp
  GIT_REPOSITORY https://github.com/jonathonl/omp.git
  GIT_TAG v1.0.0
)
FetchContent_MakeAvailable(omp)


find_package(Threads)
find_package(savvy REQUIRED)

add_executable(minimac4
        src/dosage_writer.cpp src/dosage_writer.hpp
        src/hidden_markov_model.cpp src/hidden_markov_model.hpp
        src/input_prep.cpp src/input_prep.hpp
        src/recombination.cpp src/recombination.hpp
        src/unique_haplotype.cpp src/unique_haplotype.hpp
        src/main.cpp)

target_link_libraries(minimac4 savvy Threads::Threads omp)
target_include_directories(minimac4 PRIVATE ${omp_SOURCE_DIR}/include)

add_custom_target(manuals
                  COMMAND help2man --output "${CMAKE_BINARY_DIR}/minimac4.1" "${CMAKE_BINARY_DIR}/minimac4")

if(BUILD_TESTS)
    enable_testing()
    add_test(simple_test ${CMAKE_SOURCE_DIR}/test/simple-test.sh ${CMAKE_BINARY_DIR}/minimac4 m4_simple_test_output ${CMAKE_SOURCE_DIR}/test/input/ref_panel.vcf ${CMAKE_SOURCE_DIR}/test/input/tar_panel.vcf)
    add_test(haploid_and_diploid_test ${CMAKE_SOURCE_DIR}/test/simple-test.sh ${CMAKE_BINARY_DIR}/minimac4 m4_haploid_and_diploid_test_output ${CMAKE_SOURCE_DIR}/test/input/ref_panel_with_haploid_sample.vcf ${CMAKE_SOURCE_DIR}/test/input/tar_panel_with_haploid_sample.vcf)
endif()

install(TARGETS minimac4 RUNTIME DESTINATION bin)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/minimac4.1
        DESTINATION share/man/man1
        OPTIONAL)

set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
include(CPack)
