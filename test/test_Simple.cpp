#include <gtest/gtest.h>
#include "run_main.hpp"

#ifndef TEST_DATA
#define TEST_DATA
#endif

TEST(Simple_run, compress)
{
    // Create args string
    std::vector<std::string> compress_args{
        "--compress-reference", std::string(TEST_DATA) + "/ref_panel.vcf.gz"
    };

    run_imputation_test(compress_args);
}