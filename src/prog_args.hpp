#ifndef MINIMAC4_PROG_ARGS_HPP
#define MINIMAC4_PROG_ARGS_HPP

#include "getopt_wrapper.hpp"

#include <savvy/reader.hpp>

#include <string>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iterator>

/**
 * @brief Stores and manages program arguments for Minimac4.
 *
 * This class wraps parsed command-line options and provides
 * typed accessors for downstream modules.
 */
class prog_args : public getopt_wrapper
{
private:
  std::string ref_path_;               ///< Path to the reference panel file.
  std::string tar_path_;               ///< Path to the target panel file.
  std::string map_path_;               ///< Path to the genetic map file.
  std::string out_path_ = "/dev/stdout"; ///< Path to the output file (default: stdout).
  std::string temp_prefix_;            ///< Prefix for temporary files.
  std::string prefix_;                 ///< Deprecated: old prefix option.
  std::string emp_out_path_;           ///< Path for empirical R2 output.
  std::string sites_out_path_;         ///< Path for sites-only output.
  savvy::file::format out_format_ = savvy::file::format::sav; ///< Output file format.
  std::uint8_t out_compression_ = 6;   ///< Compression level for output file.
  std::vector<std::string> fmt_fields_ = {"HDS"}; ///< FORMAT fields to include in output.
  std::unordered_set<std::string> sample_ids_; ///< Subset of sample IDs to process.
  savvy::genomic_region reg_ = {""};   ///< Genomic region to restrict processing.
  std::size_t temp_buffer_ = 200;      ///< Buffer size for temporary storage.
  std::size_t min_block_size_ = 10;    ///< Minimum block size for imputation.
  std::size_t max_block_size_ = 0xFFFF;///< Maximum block size (S1R limit).
  std::size_t slope_unit_ = 10;        ///< Unit size for slope calculation.
  std::int64_t chunk_size_ = 20000000; ///< Size of chunks to process (bp).
  std::int64_t overlap_ = 3000000;     ///< Overlap between chunks (bp).
  std::int16_t threads_ = 1;           ///< Number of computation threads.
  float decay_ = 0.f;                  ///< Decay parameter for HMM.
  float min_r2_ = -1.f;                ///< Minimum imputation R2 threshold.
  float min_ratio_ = 1e-4f;            ///< Minimum ratio for haplotype pruning.
  float prob_threshold_ = 0.01f;       ///< Probability threshold for genotype calls.
  float prob_threshold_s1_ = -1.f;     ///< Probability threshold for S1 records.
  float diff_threshold_ = 0.01f;       ///< Threshold for likelihood differences.
  float min_recom_ = 1e-5f;            ///< Minimum recombination rate.
  float error_param_ = 0.01f;          ///< Genotyping error parameter.
  bool all_typed_sites_ = false;       ///< Process all typed sites if true.
  bool update_m3vcf_ = false;          ///< Update M3VCF reference if true.
  bool compress_reference_ = false;    ///< Compress reference panel if true.
  bool pass_only_ = false;             ///< Keep only PASS variants if true.
  bool meta_ = false;                  ///< Deprecated: meta option.
  bool fail_min_ratio_ = true;         ///< Whether to fail if min ratio not met.
  bool help_ = false;                  ///< Show help and exit if true.
  bool version_ = false;               ///< Show version and exit if true.

public:
  /** @return true if `--help` flag is set. */
  bool help_is_set() const { return help_; }

  /** @return true if `--version` flag is set. */
  bool version_is_set() const { return version_; }

  /** @return Reference panel path. */
  const std::string& ref_path() const { return ref_path_; }

  /** @return Target panel path. */
  const std::string& tar_path() const { return tar_path_; }

  /** @return Genetic map path. */
  const std::string& map_path() const { return map_path_; }

  /** @return Output file path. */
  const std::string& out_path() const { return out_path_; }

  /** @return Empirical R2 output path. */
  const std::string& emp_out_path() const { return emp_out_path_; }

  /** @return Sites-only output path. */
  const std::string& sites_out_path() const { return sites_out_path_; }

  /** @return Prefix for temporary files. */
  const std::string& temp_prefix() const { return temp_prefix_; }

  /** @return Output file format. */
  savvy::file::format out_format() const { return out_format_; }

  /** @return Compression level for output. */
  std::uint8_t out_compression() const { return out_compression_; }

  /** @return FORMAT fields to output. */
  const std::vector<std::string>& fmt_fields() const { return fmt_fields_; }

  /** @return Subset of sample IDs. */
  const std::unordered_set<std::string>& sample_ids() const { return sample_ids_; }

  /** @return Genomic region to restrict analysis. */
  const savvy::genomic_region& region() const { return reg_; }

  /** @return Chunk size in base pairs. */
  std::int64_t chunk_size() const { return chunk_size_; }

  /** @return Overlap size in base pairs. */
  std::int64_t overlap() const { return overlap_; }

  /** @return Number of threads to use. */
  std::int16_t threads() const { return threads_; }

  /** @return Temporary buffer size. */
  std::size_t temp_buffer() const { return temp_buffer_ ; }

  /** @return Minimum block size. */
  std::size_t min_block_size() const { return min_block_size_; }

  /** @return Maximum block size. */
  std::size_t max_block_size() const { return max_block_size_; }

  /** @return Slope unit size. */
  std::size_t slope_unit() const { return slope_unit_; }

  /** @return HMM decay parameter. */
  float decay() const { return decay_; }

  /** @return Minimum imputation R2 threshold. */
  float min_r2() const { return min_r2_; }

  /** @return Minimum haplotype pruning ratio. */
  float min_ratio() const { return min_ratio_; }

  /** @return Probability threshold for genotype calls. */
  float prob_threshold() const { return prob_threshold_; }

  /** @return Probability threshold for S1 records. */
  float prob_threshold_s1() const { return prob_threshold_s1_; }

  /** @return Likelihood difference threshold. */
  float diff_threshold() const { return diff_threshold_; }

  /** @return Minimum recombination rate. */
  float min_recom() const { return min_recom_; }

  /** @return Genotyping error parameter. */
  float error_param() const { return error_param_; }

  /** @return true if all typed sites should be processed. */
  bool all_typed_sites() const { return all_typed_sites_; }

  /** @return true if updating M3VCF reference is enabled. */
  bool update_m3vcf() const { return update_m3vcf_; }

  /** @return true if reference panel should be compressed. */
  bool compress_reference() const { return compress_reference_; }

  /** @return true if only PASS variants are kept. */
  bool pass_only() const { return pass_only_; }

  /** @return true if failing on min ratio violation is enabled. */
  bool fail_min_ratio() const { return fail_min_ratio_; }

  /**
   * @brief Construct program arguments parser for minimac4.
   *
   * This constructor initializes the getopt-based argument parser with usage
   * instructions and all supported command-line options for minimac4.
   *
   * Usage:
   * @code
   *   minimac4 [opts ...] <reference.msav> <target.{sav,bcf,vcf.gz}>
   *   minimac4 [opts ...] --update-m3vcf <reference.m3vcf.gz>
   *   minimac4 [opts ...] --compress-reference <reference.{sav,bcf,vcf.gz}>
   * @endcode
   *
   * Supported options include:
   * - General:
   *   - `--help, -h` : Print usage.
   *   - `--version, -v` : Print version.
   * - Input/Output:
   *   - `--output, -o <path>` : Output path (default: /dev/stdout).
   *   - `--output-format, -O <fmt>` : Output format (bcf, sav, vcf.gz, …; default: sav).
   *   - `--sites, -s <path>` : Output path for sites-only file.
   *   - `--empirical-output, -e <path>` : Path for empirical dosages.
   * - Reference/Target control:
   *   - `--all-typed-sites, -a` : Include sites that exist only in target VCF.
   *   - `--map, -m <file>` : Genetic map file.
   *   - `--region, -r <region>` : Genomic region to impute.
   *   - `--sample-ids <list>` / `--sample-ids-file <file>` : Subset samples from reference panel.
   * - Performance:
   *   - `--threads, -t <int>` : Number of threads (default: 1).
   *   - `--temp-buffer, -b <int>` : Number of samples to buffer before writing (default: 200).
   *   - `--chunk, -c <bp>` : Maximum chunk length in base pairs (default: 20,000,000).
   *   - `--overlap, -w <bp>` : Size of flanking overlap (default: 3,000,000).
   * - HMM/Imputation parameters:
   *   - `--match-error <float>` : Match error probability (default: 0.01).
   *   - `--min-r2 <float>` : Minimum estimated r² for output variants.
   *   - `--min-ratio <float>` : Minimum target/reference site ratio (default: 1e-4).
   *   - `--min-ratio-behavior <skip|fail>` : Behavior if ratio not met (default: fail).
   *   - `--min-recom <float>` : Minimum recombination probability (default: 1e-5).
   *   - `--prob-threshold <float>` : Probability threshold for template selection.
   *   - `--prob-threshold-s1 <float>` : Probability threshold in original state space.
   *   - `--diff-threshold <float>` : Probability diff threshold for template selection.
   *   - `--decay <float>` : Dosage decay in flanking regions (default: 0, disabled).
   * - Reference compression / conversion:
   *   - `--update-m3vcf` : Convert M3VCF to MVCF.
   *   - `--compress-reference` : Compress VCF/BCF/SAV into MVCF.
   *   - `--min-block-size <int>` : Minimum haplotype block size (default: 10).
   *   - `--max-block-size <int>` : Maximum haplotype block size (default: 65535).
   *   - `--slope-unit <int>` : Slope parameter for compression heuristic (default: 10).
   *
   * @note
   * Some options are deprecated but retained for backward compatibility 
   * (e.g., `--allTypedSites`, `--rsid`, `--meta`, `--refHaps`, …).
   *
   * @see getopt_wrapper
   */
  prog_args() :
    getopt_wrapper(
      "Usage: minimac4 [opts ...] <reference.msav> <target.{sav,bcf,vcf.gz}>\n"
      "       minimac4 [opts ...] --update-m3vcf <reference.m3vcf.gz>\n"
      "       minimac4 [opts ...] --compress-reference <reference.{sav,bcf,vcf.gz}>",
      {
        {"all-typed-sites", no_argument, 0, 'a', "Include in the output sites that exist only in target VCF"},
        {"temp-buffer", required_argument, 0, 'b', "Number of samples to impute before writing to temporary files (default: 200)"},
        {"chunk", required_argument, 0, 'c', "Maximum chunk length in base pairs to impute at once (default: 20000000"},
        {"empirical-output", required_argument, 0, 'e', "Output path for empirical dosages"},
        {"help", no_argument, 0, 'h', "Print usage"},
        {"format", required_argument, 0, 'f', "Comma-separated list of format fields to generate (GT, HDS, DS, GP, or SD; default: HDS)"},
        {"map", required_argument, 0, 'm', "Genetic map file"},
        {"output", required_argument, 0, 'o', "Output path (default: /dev/stdout)"},
        {"output-format", required_argument, 0, 'O', "Default output file format used for ambiguous filenames (bcf, sav, vcf.gz, ubcf, usav, or vcf; default: sav)"},
        //{"pass-only", no_argument, 0, 'p', "Only imports variants with FILTER column set to PASS"},
        {"region", required_argument, 0, 'r', "Genomic region to impute"},
        {"sites", required_argument, 0, 's', "Output path for sites-only file"},
        {"threads", required_argument, 0, 't', "Number of threads (default: 1)"},
        {"version", no_argument, 0, 'v', "Print version"},
        {"overlap", required_argument, 0, 'w', "Size (in base pairs) of overlap before and after impute region to use as input to HMM (default: 3000000)"},
        {"decay", required_argument, 0, '\x02', "Decay rate for dosages in flanking regions (default: disabled with 0)"},
        {"min-r2", required_argument, 0, '\x02', "Minimum estimated r-square for output variants"},
        {"min-ratio", required_argument, 0, '\x02', "Minimum ratio of number of target sites to reference sites (default: 1e-4)"},
        {"min-ratio-behavior", required_argument, 0, '\x02', "Behavior for when --min-ratio is not met (\"skip\" or \"fail\"; default: fail)"}, // maybe add "warn"
        {"match-error", required_argument, 0, '\x02', "Error parameter for HMM match probabilities (default: 0.01)"},
        {"min-recom", required_argument, 0, '\x02', "Minimum recombination probability (default: 1e-5)"},
        {"prob-threshold", required_argument, 0, '\x02', "Probability threshold used for template selection"},
        {"prob-threshold-s1", required_argument, 0, '\x02', "Probability threshold used for template selection in original state space"},
        {"diff-threshold", required_argument, 0, '\x02', "Probability diff threshold used in template selection"},
        {"sample-ids", required_argument, 0, '\x02', "Comma-separated list of sample IDs to subset from reference panel"},
        {"sample-ids-file", required_argument, 0, '\x02', "Text file containing sample IDs to subset from reference panel (one ID per line)"},
        {"temp-prefix", required_argument, 0, '\x02', "Prefix path for temporary output files (default: ${TMPDIR}/m4_)"},
        {"update-m3vcf", no_argument, 0, '\x01', "Converts M3VCF to MVCF (default output: /dev/stdout)"},
        {"compress-reference", no_argument, 0, '\x01', "Compresses VCF to MVCF (default output: /dev/stdout)"},
        {"min-block-size", required_argument, 0, '\x02', "Minimium block size for unique haplotype compression (default: 10)"},
        {"max-block-size", required_argument, 0, '\x02', "Maximum block size for unique haplotype compression (default: 65535)"},
        {"slope-unit", required_argument, 0, '\x02', "Parameter for unique haplotype compression heuristic (default: 10)"},
        // vvvv deprecated vvvv //
        {"allTypedSites", no_argument, 0, '\x01', nullptr},
        {"rsid", no_argument, 0, '\x01', nullptr},
        //{"passOnly", no_argument, 0, '\x01', nullptr},
        {"meta", no_argument, 0, '\x01', nullptr},
        {"noPhoneHome", no_argument, 0, '\x01', nullptr},
        {"referenceEstimates", no_argument, 0, '\x01', nullptr},
        {"haps", required_argument, 0, '\x02', nullptr},
        {"refHaps", required_argument, 0, '\x02', nullptr},
        {"prefix", required_argument, 0, '\x02', nullptr},
        {"mapFile", required_argument, 0, '\x02', nullptr},
        {"chr", required_argument, 0, '\x02', nullptr},
        {"start", required_argument, 0, '\x02', nullptr},
        {"end", required_argument, 0, '\x02', nullptr},
        {"window", required_argument, 0, '\x02', nullptr},
        {"ChunkOverlapMb", required_argument, 0, '\x02', nullptr},
        {"ChunkLengthMb", required_argument, 0, '\x02', nullptr},
        {"cpus", required_argument, 0, '\x02', nullptr},
        {"minRatio", required_argument, 0, '\x02', nullptr}
      })
  {
  }

  /**
   * @brief Parse command-line arguments for minimac4.
   *
   * This function processes both short and long options (via getopt_long) and 
   * positional arguments to configure program execution. It updates internal
   * state variables based on recognized options and validates argument counts.
   *
   * Behavior:
   * - Parses @p argc and @p argv using @c getopt_long().
   * - Recognizes both short (-o) and long (--output) forms.
   * - Updates internal members (e.g., file paths, numeric thresholds, flags).
   * - Handles deprecated options with warnings and backward-compatible behavior.
   * - Validates required positional arguments (reference and target files).
   * - Automatically assigns temporary prefixes and output paths if not provided.
   *
   * Command-line structure:
   * @code
   *   minimac4 [options] <reference.msav> <target.{sav,bcf,vcf.gz}>
   *   minimac4 [options] --update-m3vcf <reference.m3vcf.gz>
   *   minimac4 [options] --compress-reference <reference.{sav,bcf,vcf.gz}>
   * @endcode
   *
   * Return conditions:
   * - Returns @c true if parsing succeeds (including early exits for --help or --version).
   * - Returns @c false if:
   *   - An invalid option value is provided,
   *   - An unknown option is encountered,
   *   - The number of positional arguments is invalid.
   *
   * Side effects:
   * - Writes warnings to @c std::cerr for deprecated options.
   * - Writes error messages to @c std::cerr on invalid usage.
   * - Sets defaults for temporary prefix (@c TMPDIR or /tmp) and output suffixes.
   * - Ensures "HDS" is included in format fields if @c --empirical-output is requested.
   *
   * @param argc Number of command-line arguments.
   * @param argv Array of argument strings.
   *
   * @return @c true if parsing and validation succeed, @c false otherwise.
   *
   * @note
   *   - Recognized options include general flags (help, version), I/O control,
   *     imputation parameters, and reference compression options.
   *   - Deprecated options are retained for backward compatibility, but emit warnings.
   *
   * @see getopt_long, prog_args
   */
  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt;
    while ((opt = getopt_long(argc, argv, short_opt_string_.c_str(), long_options_.data(), &long_index)) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'a':
        all_typed_sites_ = true;
        break;
      case 'b':
        temp_buffer_ = std::size_t(std::atoll(optarg ? optarg : ""));
        break;
      case 'c':
        chunk_size_ = std::atoll(optarg ? optarg : "");
        break;
      case 'e':
        emp_out_path_ = optarg ? optarg : "";
        break;
      case 'h':
        help_ = true;
        return true;
      case 'f':
        {
          fmt_fields_ = split_string_to_vector(optarg ? optarg : "", ',');
          std::unordered_set<std::string> allowed = {"GT", "GP", "DS", "HDS", "SD"};
          for (auto it = fmt_fields_.begin(); it != fmt_fields_.end(); ++it)
          {
            if (allowed.find(*it) == allowed.end())
              return std::cerr << "Error: Invalid --format option (" << *it << ")\n", false;
          }
          break;
        }
      case 'm':
        map_path_ = optarg ? optarg : "";
        break;
      case 'o':
        out_path_ = optarg ? optarg : "";
        break;
      case 'O':
        {
          using fmt = savvy::file::format;
          std::string ot = optarg ? optarg : "";
          if (ot == "vcf")
          {
            out_format_ = fmt::vcf;
            out_compression_ = 0;
          }
          else if (ot == "vcf.gz")
          {
            out_format_ = fmt::vcf;
          }
          else if (ot == "bcf")
          {
            out_format_ = fmt::bcf;
          }
          else if (ot == "ubcf")
          {
            out_format_ = fmt::bcf;
            out_compression_ = 0;
          }
          else if (ot == "sav")
          {
            out_format_ = fmt::sav;
          }
          else if (ot == "usav")
          {
            out_format_ = fmt::sav;
            out_compression_ = 0;
          }
          else
          {
            std::cerr << "Invalid --output-format: " << ot << std::endl;
            return false;
          }
          break;
        }
      case 'p':
        pass_only_ = true;
        break;
      case 'r':
        reg_ = string_to_region(optarg ? optarg : "");
        break;
      case 's':
        sites_out_path_ = optarg ? optarg : "";
        break;
      case 't':
        threads_ = atoi(optarg ? optarg : "");
        break;
      case 'v':
        version_ = true;
        return true;
      case 'w':
        overlap_ = std::atoll(optarg ? optarg : "");
        break;
      case '\x01':
        if (std::string(long_options_[long_index].name) == "update-m3vcf")
        {
          update_m3vcf_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "compress-reference")
        {
          compress_reference_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "allTypedSites")
        {
          std::cerr << "Warning: --allTypedSites is deprecated in favor of --all-typed-sites\n";
          all_typed_sites_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "rsid")
        {
          std::cerr << "Warning: --rsid is deprecated (on by default)\n";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "passOnly")
        {
          std::cerr << "Warning: --passOnly is deprecated in favor of --pass-only\n";
          pass_only_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "meta")
        {
          std::cerr << "Warning: --meta is deprecated in favor of --empirical-output\n";
          meta_ = true;
          break;
        }
        else if (std::string(long_options_[long_index].name) == "noPhoneHome")
        {
          std::cerr << "Warning: --noPhoneHome is deprecated and ignored\n";
          break;
        }
        else if (std::string(long_options_[long_index].name) == "referenceEstimates")
        {
          std::cerr << "Warning: --referenceEstimates is deprecated and ignored\n";
          break;
        }
        // else pass through to default
      case '\x02':
        {
          std::string long_opt_str = std::string(long_options_[long_index].name);
          if (long_opt_str == "decay")
          {
            decay_ = std::atof(optarg ? optarg : "");
            break;
          }
          else if (long_opt_str == "min-r2")
          {
            min_r2_ = std::atof(optarg ? optarg : "");
            break;
          }
          else if (long_opt_str == "min-ratio")
          {
            min_ratio_ = std::min(1., std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "min-ratio-behavior")
          {
            fail_min_ratio_ = std::string(optarg ? optarg : "") == "fail";
            break;
          }
          else if (long_opt_str == "match-error")
          {
            error_param_ = std::min(0.5, std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "min-recom")
          {
            min_recom_ = std::min(0.5, std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "prob-threshold")
          {
            prob_threshold_ = std::min(1., std::max(0., std::atof(optarg ? optarg : "")));
            break;
          }
          else if (long_opt_str == "prob-threshold-s1")
          {
            prob_threshold_s1_ = std::min(1., std::atof(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "temp-prefix")
          {
            temp_prefix_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "diff-threshold")
          {
            diff_threshold_ = std::max(0., std::atof(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "sample-ids")
          {
            auto tmp_ids = split_string_to_vector(optarg ? optarg : "", ',');
            sample_ids_.insert(tmp_ids.begin(), tmp_ids.end());
            break;
          }
          else if (long_opt_str == "sample-ids-file")
          {
            std::ifstream ifs(optarg ? optarg : "");
            sample_ids_.insert(std::istream_iterator<std::string>(ifs), std::istream_iterator<std::string>());
            break;
          }
          else if (long_opt_str == "min-block-size")
          {
            min_ratio_ = std::max(1ll, std::atoll(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "max-block-size")
          {
            min_ratio_ = std::max(1ll, std::atoll(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "slope-unit")
          {
            min_ratio_ = std::max(1ll, std::atoll(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "haps")
          {
            std::cerr << "Warning: --haps is deprecated\n";
            tar_path_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "refHaps")
          {
            std::cerr << "Warning: --refHaps is deprecated\n";
            ref_path_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "chr")
          {
            std::cerr << "Warning: --chr is deprecated in favor of --region\n";
            reg_ = savvy::genomic_region(optarg ? optarg : "", reg_.from(), reg_.to());
            break;
          }
          else if (long_opt_str == "start")
          {
            std::cerr << "Warning: --start is deprecated in favor of --region\n";
            reg_ = savvy::genomic_region(reg_.chromosome(), std::atoll(optarg ? optarg : ""), reg_.to());
            break;
          }
          else if (long_opt_str == "end")
          {
            std::cerr << "Warning: --end is deprecated in favor of --region\n";
            reg_ = savvy::genomic_region(reg_.chromosome(), reg_.from(), std::atoll(optarg ? optarg : ""));
            break;
          }
          else if (long_opt_str == "prefix")
          {
            std::cerr << "Warning: --prefix is deprecated in favor of --output, --empirical-output, and --sites\n";
            prefix_ = optarg ? optarg : "";
            out_format_ = savvy::file::format::vcf; // Default to VCF when --prefix is used in order to be consistent with previous behavior.
            out_compression_ = 6;
            break;
          }
          else if (long_opt_str == "mapFile")
          {
            std::cerr << "Warning: --mapFile is deprecated in favor of --map\n";
            map_path_ = optarg ? optarg : "";
            break;
          }
          else if (long_opt_str == "window")
          {
            std::cerr << "Warning: --window is deprecated in favor of --overlap\n";
            overlap_ = std::atoll(optarg ? optarg : "");
            break;
          }
          else if (long_opt_str == "ChunkLengthMb")
          {
            std::cerr << "Warning: --ChunkLengthMb is deprecated in favor of --chunk\n";
            chunk_size_ = std::atoll(optarg ? optarg : "") * 1000000;
            break;
          }
          else if (long_opt_str == "ChunkOverlapMb")
          {
            std::cerr << "Warning: --ChunkOverlapMb is deprecated in favor of --overlap\n";
            overlap_ = std::atoll(optarg ? optarg : "") * 1000000;
            break;
          }
          else if (long_opt_str == "cpus")
          {
            std::cerr << "Warning: --cpus is deprecated in favor of --threads\n";
            threads_ = atoi(optarg ? optarg : "");
            break;
          }
          else if (long_opt_str == "minRatio")
          {
            std::cerr << "Warning: --minRatio is deprecated in favor of --min-ratio\n";
            min_ratio_ = std::atof(optarg ? optarg : "");
            break;
          }
          // else pass through to default
        }
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 2)
    {
      ref_path_ = argv[optind];
      tar_path_ = argv[optind + 1];
    }
    else if ((update_m3vcf_ || compress_reference_) && remaining_arg_count == 1)
    {
      ref_path_ = argv[optind];
    }
    else if (remaining_arg_count < 2)
    {
      if (ref_path_.empty() || tar_path_.empty())
      {
        std::cerr << "Too few arguments\n";
        return false;
      }
    }
    else
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (!prefix_.empty())
    {
      std::string suffix = "sav";
      if (out_format_ == savvy::file::format::bcf)
        suffix = "bcf";
      else if (out_format_ == savvy::file::format::vcf)
      {
        suffix = "vcf";
        if (out_compression_)
          suffix += ".gz";
      }

      out_path_ = prefix_ + ".dose." + suffix;
      sites_out_path_ = prefix_ + ".sites." + suffix;
      if (meta_)
        emp_out_path_ = prefix_ + ".empiricalDose." + suffix;
    }

    if (temp_prefix_.empty())
    {
      char* tmpdir = std::getenv("TMPDIR");
      if (tmpdir && strlen(tmpdir))
      {
        std::string p(tmpdir);
        if (p.back() != '/')
          p += '/';
        temp_prefix_ = p + "m4_";
      }
      else
      {
        temp_prefix_ = "/tmp/m4_";
      }
    }

    if (!emp_out_path_.empty() && std::find(fmt_fields_.begin(), fmt_fields_.end(), "HDS") == fmt_fields_.end())
      fmt_fields_.emplace_back("HDS");

    return true;
  }
private:
  /**
   * @brief Convert a genomic region string into a savvy::genomic_region object.
   *
   * This function parses a region string of the form:
   *   - `"chr"` → entire chromosome (e.g., `"chr1"`)
   *   - `"chr:pos"` → single position (e.g., `"chr1:12345"`)
   *   - `"chr:start-end"` → interval (e.g., `"chr1:1000-2000"`)
   *   - `"chr:start-"` → interval with open end (treated as chr:start only)
   *
   * Parsing logic:
   * - If no colon (':') is found, the string is interpreted as a chromosome.
   * - If a colon is found but no hyphen ('-'), the string is interpreted as a single locus.
   * - If both colon and hyphen are found:
   *   - If the end coordinate is missing, only the start coordinate is used.
   *   - Otherwise, both start and end coordinates are parsed.
   *
   * @param s Input string specifying a genomic region.
   * @return A @c savvy::genomic_region corresponding to the parsed chromosome and coordinates.
   *
   * @note
   * - Coordinates are parsed as 64-bit unsigned integers via @c std::atoll.
   * - No validation is performed on chromosome naming or coordinate ordering.
   * - If parsing fails (non-numeric positions), behavior depends on @c std::atoll (undefined / 0 fallback).
   *
   * @see savvy::genomic_region
   */
  savvy::genomic_region string_to_region(const std::string& s)
  {
    const std::size_t colon_pos = s.find(':');
    if (colon_pos == std::string::npos)
    {
      return savvy::genomic_region(s);
    }
    else
    {
      std::string chr = s.substr(0, colon_pos);
      const std::size_t hyphen_pos = s.find('-', colon_pos + 1);
      if (hyphen_pos == std::string::npos)
      {
        std::string slocus = s.substr(colon_pos + 1);
        std::uint64_t ilocus = std::uint64_t(std::atoll(slocus.c_str()));
        return savvy::genomic_region(chr, ilocus, ilocus);
      }
      else
      {
        std::string sbeg = s.substr(colon_pos + 1, hyphen_pos - chr.size() - 1);
        std::string send = s.substr(hyphen_pos + 1);
        if (send.empty())
        {
          return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())));
        }
        else
        {
          return savvy::genomic_region(chr, std::uint64_t(std::atoll(sbeg.c_str())), std::uint64_t(std::atoll(send.c_str())));
        }
      }
    }
  }

  /**
   * @brief Split a C-style string into tokens based on a delimiter.
   *
   * This function takes a null-terminated C string and splits it into
   * substrings whenever the specified delimiter character is encountered.
   *
   * Example:
   * @code
   * const char* input = "apple,banana,cherry";
   * auto tokens = split_string_to_vector(input, ',');
   * // tokens = {"apple", "banana", "cherry"}
   * @endcode
   *
   * @param in   Input null-terminated string to split.
   * @param delim Delimiter character used to split the string.
   * @return A vector of substrings, each as a std::string.
   *
   * @note
   * - Consecutive delimiters will produce empty strings in the result.
   * - If the delimiter does not appear in the string, the entire input is returned as a single element.
   * - Uses @c std::find to locate delimiters efficiently.
   */
  std::vector<std::string> split_string_to_vector(const char* in, char delim)
  {
    std::vector<std::string> ret;
    const char* d = nullptr;
    std::string token;
    const char* s = in;
    const char*const e = in + strlen(in);
    while ((d = std::find(s, e,  delim)) != e)
    {
      ret.emplace_back(std::string(s, d));
      s = d ? d + 1 : d;
    }
    ret.emplace_back(std::string(s,d));
    return ret;
  }
};

#endif // MINIMAC4_PROG_ARGS_HPP
