## Install OMP
include(FetchContent)

FetchContent_Declare(
  omp
  GIT_REPOSITORY https://github.com/jonathonl/omp.git
  GIT_TAG v1.0.0
)
FetchContent_MakeAvailable(omp)

## Find local packages
find_package(Threads)
find_package(savvy REQUIRED)

## Create source library
add_library(minimac4_source dosage_writer.cpp
                            hidden_markov_model.cpp
                            input_prep.cpp
                            recombination.cpp
                            unique_haplotype.cpp
                            imputation.cpp
)

## Link to header files
target_include_directories(minimac4_source PRIVATE ${omp_SOURCE_DIR}/include)

## Link to dependencies
target_link_libraries(minimac4_source PRIVATE savvy Threads::Threads omp)