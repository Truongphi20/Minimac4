# Invoke Gtest
fetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

## Simple test
add_executable(test_Simple_compress test_Simple_compress.cpp run_main.cpp)
target_link_libraries(test_Simple_compress GTest::gtest_main minimac4_source)
target_compile_definitions(test_Simple_compress PRIVATE TEST_DATA="${CMAKE_SOURCE_DIR}/test/data")
add_test(NAME test_Simple_compress COMMAND test_Simple_compress)

add_executable(test_Simple_impute test_Simple_impute.cpp run_main.cpp)
target_link_libraries(test_Simple_impute GTest::gtest_main minimac4_source)
target_compile_definitions(test_Simple_impute PRIVATE TEST_DATA="${CMAKE_SOURCE_DIR}/test/data")
add_test(NAME test_Simple_impute COMMAND test_Simple_impute)

