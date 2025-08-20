#pragma once

#include <cstdio>
#include <iostream>
#include <ctime>
#include <unistd.h>
#include <fstream>
#include <ostream>

#include "prog_args.hpp"
#include "input_prep.hpp"
#include "hidden_markov_model.hpp"
#include "recombination.hpp"
#include "dosage_writer.hpp"

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>
#include <omp.hpp>

/**
 * @class imputation
 * @brief Class responsible for managing genotype imputation statistics and timing.
 *
 * This class tracks cumulative input, output, and imputation times for 
 * genotype imputation runs. It also provides helper methods to record 
 * elapsed times and retrieve totals.
 *
 * The timing values are updated by internal helper functions whenever 
 * an operation (e.g., input loading, output writing, or HMM traversal) 
 * completes. These values can then be retrieved for reporting and profiling.
 */
class imputation
{
    private:
        /**
     * @brief Accumulated total time spent loading input (in seconds).
     */
    long total_input_time_ = 0;

    /**
     * @brief Accumulated total time spent writing output (in seconds).
     */
    long total_output_time_ = 0;

    /**
     * @brief Accumulated total time spent on imputation (in seconds).
     */
    long total_impute_time_ = 0;
    private:
        /**
         * @brief Record elapsed input time and update cumulative total.
         * @param diff Elapsed time (in seconds).
         * @return The same elapsed time value.
         */
        double record_input_time(double diff) { total_input_time_ += diff; return diff; }

        /**
         * @brief Record elapsed output time and update cumulative total.
         * @param diff Elapsed time (in seconds).
         * @return The same elapsed time value.
         */
        double record_output_time(double diff) { total_output_time_ += diff; return diff; }

        /**
         * @brief Record elapsed imputation time and update cumulative total.
         * @param diff Elapsed time (in seconds).
         * @return The same elapsed time value.
         */
        double record_impute_time(double diff) { total_impute_time_ += diff; return diff; }
    public:
        /**
         * @brief Get the total accumulated input time.
         * @return Total input time in seconds.
         */
        long total_input_time() const  { return total_input_time_; }

        /**
         * @brief Get the total accumulated output time.
         * @return Total output time in seconds.
         */
        long total_output_time() const  { return total_output_time_; }

        /**
         * @brief Get the total accumulated imputation time.
         * @return Total imputation time in seconds.
         */
        long total_impute_time() const  { return total_impute_time_; }

        /**
         * @brief Perform genotype imputation for a given genomic region.
         *
         * This function loads target haplotypes, reference haplotypes, and executes
         * a Hidden Markov Model (HMM)–based imputation algorithm across the specified
         * region. The region is optionally extended by an overlap to reduce edge effects.  
         * Results are written to the provided output writer or via temporary files 
         * when sample sizes are large.
         *
         * Workflow:
         *  - Extend the target imputation region with overlap.
         *  - Load target haplotypes for the region.
         *  - Load reference haplotypes (typed-only + full reference panel).
         *  - Check ratio of typed to reference variants against thresholds.
         *  - Run HMM forward/backward traversal in parallel for each sample haplotype.
         *  - Write dosages to temporary files (if buffered) or directly to output.
         *  - Merge temporary files if necessary.
         *
         * @param impute_region  The genomic region to impute (chromosome, start, end).
         * @param args           Program arguments controlling paths, thresholds, buffers, and options.
         * @param tpool          Thread pool for parallel HMM traversal.
         * @param output         Dosage writer for writing final imputation results.
         *
         * @return True if the imputation completed successfully (even if the chunk was skipped),
         *         false if an error occurred (e.g., file loading/writing failed).
         *
         * @note
         *  - If no reference variants exist in the region, the chunk is skipped.
         *  - If the ratio of typed sites to reference sites is below `args.min_ratio()`,
         *    the chunk may be skipped or treated as an error depending on `args.fail_min_ratio()`.
         *  - Target-only variants can be optionally included in the output (`--all-typed-sites`).
         *  - Temporary files are created and merged automatically when processing in buffered groups.
         */
        bool impute_chunk(const savvy::region& impute_region, const prog_args& args, omp::internal::thread_pool2& tpool, dosage_writer& output);
};