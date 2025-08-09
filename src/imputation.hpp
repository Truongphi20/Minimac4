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


class imputation
{
    private:
        long total_input_time_ = 0;
        long total_output_time_ = 0;
        long total_impute_time_ = 0;
    private:
        double record_input_time(double diff) { total_input_time_ += diff; return diff; }
        double record_output_time(double diff) { total_output_time_ += diff; return diff; }
        double record_impute_time(double diff) { total_impute_time_ += diff; return diff; }
    public:
        long total_input_time() const  { return total_input_time_; }
        long total_output_time() const  { return total_output_time_; }
        long total_impute_time() const  { return total_impute_time_; }

        bool impute_chunk(const savvy::region& impute_region, const prog_args& args, omp::internal::thread_pool2& tpool, dosage_writer& output);
};