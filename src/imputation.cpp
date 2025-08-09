#include "imputation.hpp"

bool imputation::impute_chunk(const savvy::region& impute_region, const prog_args& args, omp::internal::thread_pool2& tpool, dosage_writer& output)
{
    savvy::region extended_region =
        {
        impute_region.chromosome(),
        std::uint64_t(std::max(std::int64_t(1), std::int64_t(impute_region.from()) - args.overlap())),
        impute_region.to() + args.overlap()
        };

    std::cerr << "Imputing " << impute_region.chromosome() << ":" << impute_region.from() << "-" << impute_region.to() << " ..." << std::endl;

    std::cerr << "Loading target haplotypes ..." << std::endl;
    std::time_t start_time = std::time(nullptr);
    std::vector<std::string> sample_ids;
    std::vector<target_variant> target_sites;
    if (!load_target_haplotypes(args.tar_path(), extended_region, target_sites, sample_ids))
        return std::cerr << "Error: failed loading target haplotypes\n", false;
    std::cerr << "Loading target haplotypes took " << record_input_time(std::difftime(std::time(nullptr), start_time)) << " seconds" << std::endl;

    std::cerr << "Loading reference haplotypes ..." << std::endl;
    start_time = std::time(nullptr);
    reduced_haplotypes typed_only_reference_data(16, 512);
    reduced_haplotypes full_reference_data;
    std::unique_ptr<genetic_map_file> mf(args.map_path().empty() ? nullptr : new genetic_map_file(args.map_path(), impute_region.chromosome()));
    if (!load_reference_haplotypes(args.ref_path(), extended_region, impute_region, args.sample_ids(), target_sites, typed_only_reference_data, full_reference_data, mf.get(), args.min_recom(), args.error_param()))
        return std::cerr << "Error: failed loading reference haplotypes\n", false;
    std::cerr << "Loading reference haplotypes took " << record_input_time(std::difftime(std::time(nullptr), start_time)) << " seconds" << std::endl;


    std::vector<target_variant> target_only_sites = separate_target_only_variants(target_sites);

    double impute_time = 0.;
    double temp_write_time = 0.;

    std::list<savvy::reader> temp_files;
    std::list<savvy::reader> temp_emp_files;
    //    std::list<std::string> temp_files;
    //    std::list<std::string> temp_emp_files;
    full_dosages_results hmm_results;

    if (full_reference_data.variant_size() == 0)
    {
        std::cerr << "Notice: skipping empty region in reference (" << impute_region.chromosome() << ":" << impute_region.from() << ":" << impute_region.to() << ")" << std::endl;
    }
    else
    {
        float tar_ref_ratio = float(typed_only_reference_data.variant_size()) / float(full_reference_data.variant_size());
        std::cerr << "Typed sites to imputed sites ratio: " << tar_ref_ratio << " (" << typed_only_reference_data.variant_size() << "/" << full_reference_data.variant_size() << ")\n";
        if (tar_ref_ratio < args.min_ratio())
        {
        std::cerr << (args.fail_min_ratio() ? "Error" : "Warning")  << ": not enough target variants are available to impute this chunk. The --min-ratio, --chunk, or --region options may need to be altered." << std::endl;
        if (!args.fail_min_ratio())
            std::cerr << "Warning: skipping chunk " << impute_region.chromosome() << ":" << impute_region.from() << "-" << impute_region.to() << std::endl;
        if (args.fail_min_ratio())
            return false;
        return true; // skip
        }

        if (target_only_sites.size())
        {
        std::size_t cnt = std::count_if(target_only_sites.begin(), target_only_sites.end(), [impute_region](target_variant& v){ return v.pos >= impute_region.from() && v.pos <= impute_region.to(); });
        std::cerr << cnt << " variants are exclusive to target file and will be ";
        if (args.all_typed_sites())
            std::cerr << "included in output\n";
        else
        {
            std::cerr << "excluded from output\n";
            target_only_sites.clear();
        }
        }

        if (target_sites.empty())
        return std::cerr << "Error: no target variants\n", false;

    //      std::cerr << "Loading switch probabilities ..." << std::endl;
    //      start_time = std::time(nullptr);
    //      if (!load_variant_hmm_params(target_sites, typed_only_reference_data, args.error_param(), args.min_recom(), args.map_path()))
    //        return std::cerr << "Error: parsing map file failed\n", false;
    //      std::cerr << "Loading switch probabilities took " << record_input_time(std::difftime(std::time(nullptr), start_time)) << " seconds" << std::endl;

        auto reverse_maps = generate_reverse_maps(typed_only_reference_data);

        std::cerr << "Running HMM with " << args.threads() << " threads ..." << std::endl;
        start_time = std::time(nullptr);
        std::vector<hidden_markov_model> hmms(tpool.thread_count(), {args.prob_threshold(), args.prob_threshold_s1(), args.diff_threshold(), 1e-5f, args.decay()});

        std::size_t ploidy = target_sites[0].gt.size() / sample_ids.size();
        std::size_t haplotype_buffer_size = args.temp_buffer() * ploidy;
        assert(ploidy && target_sites[0].gt.size() % sample_ids.size() == 0);

        hmm_results.resize(full_reference_data.variant_size(), target_sites.size(), std::min(haplotype_buffer_size, target_sites[0].gt.size()));

        for (std::size_t i = 0; i < target_sites[0].gt.size(); i += haplotype_buffer_size)
        {
        std::size_t group_size = std::min(target_sites[0].gt.size() - i, haplotype_buffer_size);
        if (group_size < haplotype_buffer_size)
        {
            for (std::size_t j = 0; j < hmm_results.dosages_.size(); ++j)
            hmm_results.dosages_[j].resize(group_size);
            for (std::size_t j = 0; j < hmm_results.loo_dosages_.size(); ++j)
            hmm_results.loo_dosages_[j].resize(group_size);
        }

        if (i > 0)
            hmm_results.fill_eov();

        start_time = std::time(nullptr);
        omp::parallel_for_exp(
            omp::static_schedule(), omp::sequence_iterator(i), omp::sequence_iterator(i + group_size), [&](int& i, const omp::iteration_context& ctx)
            {
            if (savvy::typed_value::is_end_of_vector(target_sites[0].gt[i]))
                return; // Sample has fewer haplotypes
            hmms[ctx.thread_index].traverse_forward(typed_only_reference_data.blocks(), target_sites, i);
            hmms[ctx.thread_index].traverse_backward(typed_only_reference_data.blocks(), target_sites, i, i % haplotype_buffer_size, reverse_maps, hmm_results, full_reference_data);
            },
            tpool);
        impute_time += std::difftime(std::time(nullptr), start_time);

        start_time = std::time(nullptr);

        int tmp_fd = -1;
        int tmp_emp_fd = -1;

        if (target_sites[0].gt.size() > haplotype_buffer_size)
        {
            std::string out_emp_path;
            std::string out_path = args.temp_prefix() + std::to_string(i / haplotype_buffer_size) + "_XXXXXX";
            tmp_fd = mkstemp(&out_path[0]);
            if (tmp_fd < 0)
            return std::cerr << "Error: could not open temp file (" << out_path << ")" << std::endl, false;

            if (args.emp_out_path().size())
            {
            out_emp_path =  args.temp_prefix() + std::to_string(i / haplotype_buffer_size) + "_emp_XXXXXX";
            tmp_emp_fd = mkstemp(&out_emp_path[0]);
            if (tmp_emp_fd < 0)
                return std::cerr << "Error: could not open temp file (" << out_emp_path << ")" << std::endl, false;
            }

            dosage_writer temp_output(out_path, out_emp_path,
            "", // sites path
            savvy::file::format::sav,
            std::min<std::uint8_t>(3, args.out_compression()),
            {sample_ids.begin() + (i / ploidy), sample_ids.begin() + (i + group_size) / ploidy},
            {"HDS"},
            impute_region.chromosome(),
            -1.f, true);

            temp_files.emplace_back(out_path);
            ::close(tmp_fd);
            std::remove(out_path.c_str());
            assert(tmp_fd > 0);

            if (out_emp_path.size())
            {
            temp_emp_files.emplace_back(out_emp_path);
            ::close(tmp_emp_fd);
            std::remove(out_emp_path.c_str());
            assert(tmp_emp_fd > 0);
            }

            if (!temp_output.write_dosages(hmm_results, target_sites, target_only_sites, {i, i + group_size}, full_reference_data, impute_region))
            return std::cerr << "Error: failed writing output\n", false;
            temp_write_time += std::difftime(std::time(nullptr), start_time);

            std::cerr << "Completed " << (i + group_size) / ploidy << " of " << sample_ids.size() << " samples" << std::endl;
        }
        }

        std::cerr << "Running HMM took " << record_impute_time(impute_time) << " seconds" << std::endl;
    }

    if (temp_files.size())
    {
        std::cerr << "Writing temp files took " << record_output_time(temp_write_time) << " seconds" << std::endl;

        std::cerr << "Merging temp files ... " << std::endl;
        start_time = std::time(nullptr);
        //dosage_writer output(args.out_path(), args.emp_out_path(), args.sites_out_path(), args.out_format(), args.out_compression(), sample_ids, args.fmt_fields(), target_sites.front().chrom, false);
        if (!output.merge_temp_files(temp_files, temp_emp_files))
        return std::cerr << "Error: failed merging temp files\n", false;
        std::cerr << "Merging temp files took " << record_output_time(std::difftime(std::time(nullptr), start_time)) << " seconds" << std::endl;
    }
    else
    {
        std::size_t n_tar_haps = 0;
        if (target_sites.size())
        n_tar_haps = target_sites[0].gt.size();
        else if (target_only_sites.size())
        n_tar_haps = target_only_sites[0].gt.size(); // For when reference panel has no variants in region.

        if (n_tar_haps)
        {
        std::cerr << "Writing output ... " << std::endl;
        start_time = std::time(nullptr);
        if (!output.write_dosages(hmm_results, target_sites, target_only_sites, {0, n_tar_haps}, full_reference_data, impute_region))
            return std::cerr << "Error: failed writing output\n", false;
        std::cerr << "Writing output took " << record_output_time(std::difftime(std::time(nullptr), start_time)) << " seconds" << std::endl;
        }
    }

    std::cerr << std::endl;

    return true;
}