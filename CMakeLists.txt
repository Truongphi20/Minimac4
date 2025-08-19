cmake_minimum_required(VERSION 3.14)
project(minimac4 VERSION 4.1.6)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type (Debug|Release|RelWithDebInfo|MinSizeRel)" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)

add_definitions(-DVERSION="${PROJECT_VERSION}")

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

## Setup executor
add_subdirectory(src)

option(BUILD_TESTS "Build test cases" OFF)
if(BUILD_TESTS)
    set(CMAKE_BUILD_TYPE Debug)
    enable_testing()
    add_subdirectory(test)
    # add_test(simple_test ${CMAKE_SOURCE_DIR}/test/simple-test.sh ${CMAKE_BINARY_DIR}/minimac4 m4_simple_test_output ${CMAKE_SOURCE_DIR}/test/input/ref_panel.vcf ${CMAKE_SOURCE_DIR}/test/input/tar_panel.vcf)
    # add_test(haploid_and_diploid_test ${CMAKE_SOURCE_DIR}/test/simple-test.sh ${CMAKE_BINARY_DIR}/minimac4 m4_haploid_and_diploid_test_output ${CMAKE_SOURCE_DIR}/test/input/ref_panel_with_haploid_sample.vcf ${CMAKE_SOURCE_DIR}/test/input/tar_panel_with_haploid_sample.vcf)
else()
  add_executable(minimac4 src/main.cpp)
  target_link_libraries(minimac4 minimac4_source)
  add_custom_target(manuals
                    COMMAND help2man --output "${CMAKE_BINARY_DIR}/minimac4.1" "${CMAKE_BINARY_DIR}/minimac4")

  install(TARGETS minimac4 RUNTIME DESTINATION bin)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/minimac4.1
          DESTINATION share/man/man1
          OPTIONAL)

  set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
  set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
  set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
  include(CPack)
endif()

## Create Doxygen (option)
option(DOXY "CREATE DOXYGEN DOCUMENT" OFF)
if(DOXY)
    set(DOXYGEN_EXECUTABLE "/opt/conda/envs/minimac4-dev/bin/doxygen" CACHE FILEPATH "Path to Doxygen executable" FORCE)
    find_package(Doxygen REQUIRED)

    ## Format css
    include(FetchContent)
    FetchContent_Declare(
        doxygen-awesome-css
        URL https://github.com/jothepro/doxygen-awesome-css/archive/refs/heads/main.zip
    )
    FetchContent_MakeAvailable(doxygen-awesome-css)
    FetchContent_GetProperties(doxygen-awesome-css SOURCE_DIR AWESOME_CSS_DIR)

    # Paths for input/output
    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # Replace variables in template
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    # Create custom target
    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating Doxygen documentation"
        VERBATIM
    )
endif()