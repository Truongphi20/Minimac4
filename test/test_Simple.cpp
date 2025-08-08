#include <gtest/gtest.h>
#include "run_main.hpp"
#include <cstdio>

#ifndef TEST_DATA
#define TEST_DATA
#endif

TEST(Simple_run, compress)
{
    // Open file and redirect stdout
    FILE* fp = freopen("output.msav", "w", stdout);

    // Create args string
    std::vector<std::string> compress_args{
        "minimac4",
        "--compress-reference", std::string(TEST_DATA) + "/ref_panel.vcf.gz"
    };

    run_imputation_test(compress_args);

    fflush(stdout); // Ensure all output is flushed
    fclose(fp);     // Close the file

    // Restore stdout back to console
    freopen("/dev/tty", "w", stdout);
}