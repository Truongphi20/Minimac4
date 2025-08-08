#include "run_main.hpp"

// Helper function to run the main imputation pipeline
int run_imputation_test(std::vector<std::string> compress_args)
{
    // Convert to cstring
    std::vector<char*> cstrings;
    cstrings.reserve(compress_args.size() + 1);

    for (auto &s : compress_args) {
        char* copy = new char[s.size() + 1];
        std::strcpy(copy, s.c_str());
        cstrings.push_back(copy);
    }
    cstrings.push_back(nullptr);

    // Run main
    std::time_t start_time = std::time(nullptr);
    prog_args args;
    args.parse(static_cast<int>(compress_args.size() + 1), cstrings.data());

    if (args.help_is_set())
    {
        args.print_usage(std::cout);
        return EXIT_SUCCESS;
    }

    if (args.version_is_set())
    {
        std::cout << "minimac v" << VERSION << std::endl;
        return EXIT_SUCCESS;
    }

    std::cerr << "minimac v" << VERSION << "\n\n";

    if (args.update_m3vcf())
        return convert_old_m3vcf(args.ref_path(), args.out_path(), args.map_path())
                   ? EXIT_SUCCESS : EXIT_FAILURE;

    if (args.compress_reference())
        return compress_reference_panel(args.ref_path(), args.out_path(),
                                        args.min_block_size(), args.max_block_size(),
                                        args.slope_unit(), args.map_path())
                   ? EXIT_SUCCESS : EXIT_FAILURE;

    std::uint64_t end_pos = args.region().to();
    std::string chrom = args.region().chromosome();
    if (!stat_ref_panel(args.ref_path(), chrom, end_pos))
        return std::cerr << "Error: could not stat reference file\n", EXIT_FAILURE;

    std::vector<std::string> sample_ids;
    if (!stat_tar_panel(args.tar_path(), sample_ids))
        return std::cerr << "Error: could not stat target file\n", EXIT_FAILURE;

    dosage_writer output(args.out_path(),
                         args.emp_out_path(),
                         args.sites_out_path(),
                         args.out_format(),
                         args.out_compression(),
                         sample_ids,
                         args.fmt_fields(),
                         chrom,
                         args.min_r2(), false);

    omp::internal::thread_pool2 tpool(args.threads());
    imputation imputer;
    for (std::uint64_t chunk_start_pos = std::max<std::uint64_t>(1, args.region().from());
         chunk_start_pos <= end_pos;
         chunk_start_pos += args.chunk_size())
    {
        std::uint64_t chunk_end_pos = std::min(end_pos, chunk_start_pos + args.chunk_size() - 1ul);
        savvy::region impute_region = {chrom, chunk_start_pos, chunk_end_pos};

        if (!imputer.impute_chunk(impute_region, args, tpool, output))
            return EXIT_FAILURE;
    }

    auto total_time = long(std::difftime(std::time(nullptr), start_time));

    output.print_mean_er2(std::cerr);
    std::cerr << std::endl;
    std::fprintf(stderr, "Total time for parsing input: %ld seconds\n", imputer.total_input_time());
    std::fprintf(stderr, "Total time for HMM: %ld seconds\n", imputer.total_impute_time());
    std::fprintf(stderr, "Total time for writing output: %ld seconds\n", imputer.total_output_time());
    std::fprintf(stderr, "Total wall time (h:mm:ss): %ld:%02ld:%02ld\n",
                 total_time / 3600, (total_time % 3600) / 60, total_time % 60);

    return EXIT_SUCCESS;
}