#include <gtest/gtest.h>
#include "run_main.hpp"
#include <cstdio>

#ifndef TEST_DATA
#define TEST_DATA
#endif

TEST(Simple_run, compress)
{
    // Open file and redirect stdout
    FILE* fp = freopen("ref_panel.msav", "w", stdout);

    // Create args string
    std::vector<std::string> compress_args{
        "minimac4",
        "--compress-reference", std::string(TEST_DATA) + "/ref_panel.vcf.gz"
    };

    // Run minimac4
    run_imputation_test(compress_args);

    // Clean up
    fflush(stdout); // Ensure all output is flushed
    fclose(fp);     // Close the file

    // Restore stdout back to console
    freopen("/dev/tty", "w", stdout);
}

TEST(Simple_run, impute)
{
    // Open file and redirect stdout
    FILE* fp = freopen("imputed.vcf.gz", "w", stdout);

    // Create args string
    std::vector<std::string> impute_args{
        "minimac4",
        std::string(TEST_DATA) + "/ref_panel.msav",
        std::string(TEST_DATA) + "/tar_panel.vcf.gz",
        "-f", "GT",
        "-O", "vcf.gz",
        "--temp-buffer", "2"
    };

    // Run minimac4
    run_imputation_test(impute_args);

    // Clean up
    fflush(stdout); // Ensure all output is flushed
    fclose(fp);     // Close the file

    // Restore stdout back to console
    freopen("/dev/tty", "w", stdout);

}