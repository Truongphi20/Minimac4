# Invoke Gtest
fetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

## Simple test
add_executable(test_Simple test_Simple.cpp)
target_link_libraries(test_Simple GTest::gtest_main minimac4_source)
add_test(NAME test_Simple COMMAND test_Simple)

