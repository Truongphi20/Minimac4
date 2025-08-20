#ifndef MINIMAC4_INPUT_PREP_HPP
#define MINIMAC4_INPUT_PREP_HPP

#include "unique_haplotype.hpp"

#include <savvy/reader.hpp>

/**
 * @brief Extract sample IDs from a target panel file.
 *
 * Opens a target haplotype/genotype file using `savvy::reader`
 * and retrieves the list of sample IDs contained in the file header.
 *
 * @param tar_file_path Path to the target panel file.
 * @param sample_ids    Reference to a vector that will be filled with
 *                      the sample IDs read from the file.
 *
 * @return True if the file was successfully opened and sample IDs were
 *         retrieved, false if the file could not be opened.
 *
 * @note
 * - The function prints an error message to `stderr` if the file
 *   cannot be opened.
 * - The existing contents of `sample_ids` are overwritten.
 */
bool stat_tar_panel(const std::string& tar_file_path, std::vector<std::string>& sample_ids);

/**
 * @brief Inspect a reference panel file to determine chromosome and end position.
 *
 * This function checks the index files associated with a reference panel 
 * (S1R, CSI, or TBI) to determine the contig (chromosome) and maximum position.  
 * It supports both `.s1r` index statistics and VCF/BCF headers when CSI/TBI indexes are present.
 *
 * @param ref_file_path Path to the reference panel file (VCF/BCF/MVCF).
 * @param chrom         Chromosome name. If empty, it will be set automatically.  
 *                      If non-empty, the function verifies that the reference file contains it. 
 * @param end_pos       End position of the region. Updated to the minimum of its current 
 *                      value and the chromosome length / max position found in the index/header.
 *
 * @return True if chromosome and end position were successfully determined,  
 *         false if an error occurred (e.g., multiple chromosomes present, 
 *         index missing, or inconsistent contig name).
 *
 * @note
 * - If `.s1r` index statistics are available, they take priority.
 * - If multiple contigs are present and `chrom` is empty, the function fails 
 *   and suggests using `--region`.
 * - If `.csi` or `.tbi` index is found, chromosome information is read 
 *   from the VCF/BCF header.
 * - For MVCF files, the reference must be indexed. Legacy M3VCF files need 
 *   conversion with:
 *   @code
 *   minimac4 --update-m3vcf input.m3vcf.gz > output.msav
 *   @endcode
 *
 * @warning
 * Prints descriptive error messages to `stderr` if reference panel validation fails.
 *
 * @dot
 * digraph stat_ref_panel_flow {
 *   node [shape=box, style=rounded];
 *
 *   start [label="Start: ref_file_path"];
 *   check_s1r [label="Check .s1r index exists?"];
 *   s1r_found [label="Use s1r::stat_index()"];
 *   chrom_set [label="chrom empty?"];
 *   chrom_ok [label="Set chrom + end_pos"];
 *   chrom_verify [label="Verify contig matches chrom"];
 *   multi_chrom [label="Multiple contigs → require --region", shape=diamond];
 *
 *   check_csi [label="Check .csi or .tbi index?"];
 *   parse_header [label="Open file, read header"];
 *   chrom_empty [label="chrom empty?"];
 *   parse_first_var [label="Read first variant to set chrom"];
 *   parse_contig [label="Parse contig length from header"];
 *   fail_parse [label="Fail: cannot parse contig length"];
 *
 *   fail [label="Fail: could not load reference index"];
 *   success [label="Success", shape=oval, style=filled, fillcolor=lightgreen];
 *
 *   start -> check_s1r;
 *   check_s1r -> s1r_found [label="yes"];
 *   check_s1r -> check_csi [label="no"];
 *   s1r_found -> chrom_set;
 *   chrom_set -> chrom_ok [label="yes"];
 *   chrom_set -> chrom_verify [label="no"];
 *   chrom_verify -> success [label="contig matches"];
 *   chrom_verify -> fail [label="contig mismatch"];
 *   chrom_ok -> success;
 *   s1r_found -> multi_chrom [label="multiple contigs"];
 *   multi_chrom -> fail;
 *
 *   check_csi -> parse_header [label="yes"];
 *   check_csi -> fail [label="no"];
 *   parse_header -> chrom_empty;
 *   chrom_empty -> parse_first_var [label="yes"];
 *   chrom_empty -> parse_contig [label="no"];
 *   parse_first_var -> parse_contig;
 *   parse_contig -> success [label="parsed OK"];
 *   parse_contig -> fail_parse [label="no length"];
 *   fail_parse -> fail;
 * }
 * @enddot
 */
bool stat_ref_panel(const std::string& ref_file_path, std::string& chrom, std::uint64_t& end_pos);

/**
 * @brief Load haplotypes from a target file for a given genomic region.
 *
 * This function opens a VCF/BCF target file, extracts sample IDs, enforces ploidy 
 * consistency across all variants, and fills the list of target variants (`target_sites`) 
 * with genotypes encoded per allele.
 *
 * @param file_path    Path to the target VCF/BCF file.  
 *                     Must be bgzipped and indexed (CSI/TBI).
 * @param reg          Genomic region to query (chromosome, start, end).  
 *                     Bounds are applied using `savvy::reader::reset_bounds()`.  
 *                     If querying fails, the function returns false.
 * @param target_sites Output vector that will be filled with `target_variant` objects, 
 *                     each representing one ALT allele at a site.  
 *                     For multi-allelic sites, one entry per ALT allele is created.
 * @param sample_ids   Output vector of sample IDs extracted from the target file header.
 *
 * @return True if the haplotypes were successfully loaded and ploidy checks passed,  
 *         false otherwise. On failure, descriptive error messages are printed to `stderr`.
 *
 * @details
 * - The function extracts genotypes (`GT` field) and converts them into binary encoding:
 *   - **Biallelic sites**: GT values are stored directly in `target_variant::gt`.
 *   - **Multiallelic sites**: For each ALT allele, a separate `target_variant` is created.  
 *     The genotype array is recoded to indicate presence (1) or absence (0) of the allele.
 * - Ploidy consistency is checked:
 *   - On first variant, sample ploidies are initialized.
 *   - On subsequent variants, `check_ploidies()` ensures all samples retain the same ploidy.
 *   - If a sample’s ploidy changes mid-file, the function errors out.
 *   - Special warning is printed for chromosome X (`X` or `chrX`) due to PAR/non-PAR handling.
 * - For missing genotype values, `quiet_NaN()` is used as a placeholder for dosage fields.
 *
 * @note
 * - Requires that the input file is properly indexed. Otherwise, region queries will fail.
 * - If the function encounters ploidy changes, the run is aborted to avoid downstream 
 *   phasing/imputation errors.
 * - End users should split chromosome X into **PAR** and **non-PAR** regions before imputation.
 *
 * @warning
 * Error messages are printed directly to `stderr`.  
 * This includes:
 * - File access errors (cannot open target file).  
 * - Region query errors (file not indexed or invalid `reg`).  
 * - Ploidy inconsistency across samples.  
 */
bool load_target_haplotypes(const std::string& file_path, const savvy::genomic_region& reg, std::vector<target_variant>& target_sites, std::vector<std::string>& sample_ids);

/**
 * @brief Load and process reference haplotypes from an MVCF file.
 *
 * This function loads reference haplotypes within an extended genomic region 
 * from an MVCF (M3VCFv3/MVCFv3) file. It extracts haplotype blocks, aligns 
 * them with the provided target variants, computes allele frequencies, 
 * recombination probabilities, and compresses the resulting haplotype data 
 * into reduced representations for imputation.
 *
 * @param file_path Path to the MVCF reference haplotype file.
 * @param extended_reg Genomic region to query from the reference file (with buffer for recombination).
 * @param impute_reg Genomic region of interest for imputation (subset of @p extended_reg).
 * @param subset_ids Subset of sample IDs to extract from the reference file. 
 *                   If empty, all samples are included.
 * @param target_sites Vector of target variants to be aligned and updated with reference information 
 *                     (allele frequency, error rate, recombination probability, etc.).
 * @param typed_only_reference_data Output container for typed-only reference haplotypes, compressed.
 * @param full_reference_data Output container for full reference haplotype data across the impute region.
 * @param map_file Optional genetic map file for interpolation of centimorgan positions. 
 *                 If provided, recombination probabilities are computed from map distances.
 * @param min_recom Minimum recombination probability to enforce between adjacent variants.
 * @param default_match_error Default genotype matching error rate used when missing in the reference file.
 *
 * @return true if the haplotypes were successfully loaded and processed, 
 *         false if an error occurred (e.g., file not found, wrong format, no overlapping samples).
 *
 * @note
 * - The reference file must be indexed MVCFv3.0/M3VCFv3.0 format.
 * - If no overlapping subset samples are found, the function fails.
 * - Variants outside the imputation region are trimmed, but recombination 
 *   probabilities are still updated based on @p extended_reg.
 * - For chromosome X, ensure PAR and non-PAR regions are imputed separately.
 *
 * @see reduced_haplotypes, target_variant, savvy::reader
 */
bool load_reference_haplotypes(const std::string& file_path,
  const savvy::genomic_region& extended_reg,
  const savvy::genomic_region& impute_reg,
  const std::unordered_set<std::string>& subset_ids,
  std::vector<target_variant>& target_sites,
  reduced_haplotypes& typed_only_reference_data,
  reduced_haplotypes& full_reference_data,
  genetic_map_file* map_file,
  float min_recom,
  float default_match_error);

/**
 * @brief Loads reference haplotypes using an older recombination-based approach.
 *
 * This function reads reference haplotypes from an MVCF/M3VCF file and integrates 
 * them with a set of target variants. It populates two reduced haplotype structures:
 * one containing only variants overlapping with the target sites 
 * (`typed_only_reference_data`), and one containing all reference haplotypes 
 * within the imputation region (`full_reference_data`).
 *
 * Recombination probabilities between consecutive variants are estimated 
 * using either the centimorgan positions from a provided genetic map 
 * (`map_file`) or the reference file annotations. Allele frequencies 
 * for overlapping target variants are updated based on reference genotypes.
 *
 * @param file_path Path to the reference haplotype file (MVCF/M3VCF format).
 * @param extended_reg Genomic region specifying the extended window to load 
 *        haplotypes from (includes buffer around imputation region).
 * @param impute_reg Genomic region specifying the imputation window 
 *        (used to trim full reference haplotype blocks).
 * @param subset_ids Optional subset of sample IDs to restrict reference samples. 
 *        If empty, all samples are used.
 * @param target_sites Vector of target variants. This vector is updated with 
 *        allele frequencies, reference overlap flags, and recombination estimates.
 * @param typed_only_reference_data Output reduced haplotypes structure containing 
 *        only variants overlapping with target sites.
 * @param full_reference_data Output reduced haplotypes structure containing 
 *        all haplotype blocks overlapping the imputation region.
 * @param map_file Optional pointer to a genetic map file. If provided, used to 
 *        interpolate centimorgan distances for recombination probability calculation.
 *
 * @return True if the reference haplotypes were successfully loaded and processed, 
 *         false otherwise.
 *
 * @note The function supports both newer MVCFv3 and older M3VCF file formats. 
 *       When an invalid or incompatible file is provided, an error message 
 *       is printed and the function returns false.
 * 
 * @warning The function assumes that the reference file is properly indexed 
 *          (MVCF) or formatted (M3VCF). Errors in input files will terminate early.
 */  
bool load_reference_haplotypes_old_recom_approach(const std::string& file_path,
  const savvy::genomic_region& extended_reg,
  const savvy::genomic_region& impute_reg,
  const std::unordered_set<std::string>& subset_ids,
  std::vector<target_variant>& target_sites,
  reduced_haplotypes& typed_only_reference_data,
  reduced_haplotypes& full_reference_data,
  genetic_map_file* map_file);

/**
 * @brief Separates target-only variants from those found in the reference panel.
 *
 * This function partitions the input vector of target variants (`target_sites`) into two groups:
 *  - Variants that are present in the reference (remain in `target_sites`).
 *  - Variants that are absent from the reference (`in_ref == false`) and are moved into a 
 *    new vector (`target_only_sites`).
 *
 * The function maintains efficient memory usage by swapping elements in place rather than 
 * performing deep copies. After execution, `target_sites` will only contain variants found 
 * in the reference panel, and `target_only_sites` will contain all others.
 *
 * @param[in,out] target_sites Vector of target variants to be partitioned. 
 *                             After the call, it will contain only reference-matching variants.
 *
 * @return std::vector<target_variant>  
 *         A vector containing all target-only variants (not present in the reference).
 *
 * @note This function modifies the input vector `target_sites` by resizing it.
 * @note Order of elements may change due to the use of `std::swap`.
 *
 * @complexity O(N), where N = number of target variants.
 */  
std::vector<target_variant> separate_target_only_variants(std::vector<target_variant>& target_sites);

/**
 * @brief Loads Hidden Markov Model (HMM) parameters for target variants.
 *
 * This function initializes error rates and recombination probabilities
 * for a sequence of target variants by combining:
 *  - Default error parameters.
 *  - Per-variant error estimates (if available).
 *  - Genetic map–based recombination rates (if a genetic map file is provided).
 *
 * Specifically:
 * - Each variant’s error parameter (`err`) is set either from the reference data
 *   or from `default_error_param` if missing (NaN).
 * - Recombination probabilities (`recom`) are computed between consecutive variants
 *   using genetic map centiMorgan (cM) positions, converted to switch probabilities.
 *   If no map file is provided, previously loaded map positions in the reference
 *   data are used.
 * - The last variant is always assigned `recom = 0.0f`, ensuring no recombination
 *   at the backward traversal boundary.
 *
 * @param[in,out] tar_variants Vector of target variants whose HMM parameters (`err`, `recom`) will be filled.  
 *                             Must have the same size as the reference haplotype set.
 * @param[in,out] typed_only_reference_data Reduced reference haplotype data aligned to `tar_variants`.  
 *                                          Provides genetic map positions and optional error rates.
 * @param[in] default_error_param Default error parameter assigned if no error rate is provided.
 * @param[in] recom_min Minimum recombination probability allowed between adjacent variants.
 * @param[in] map_file_path Path to the genetic map file. If empty, recombination rates are
 *                          computed from existing positions in `typed_only_reference_data`.
 *
 * @return `true` if parameters were successfully loaded;  
 *         `false` if the variant list is empty or if the genetic map file cannot be opened.
 *
 * @note The function asserts that the number of target variants matches
 *       the number of reference variants.
 * @note Input vectors are modified in place.
 *
 * @complexity O(N), where N = number of variants.
 */
bool load_variant_hmm_params(std::vector<target_variant>& tar_variants, reduced_haplotypes& typed_only_reference_data, float default_error_param, float recom_min, const std::string& map_file_path);

/**
 * @brief Generates reverse mapping tables for reduced haplotype blocks.
 *
 * This function constructs a three-level nested vector structure
 * (`reverse_maps`) that provides, for each block and each allele state
 * in that block, the list of haplotype indices that map to it.
 *
 * The mapping is inverted from the `unique_map()` representation in each
 * block of the `typed_only_reference_data`.
 *
 * Structure of the returned vector:
 * - `reverse_maps[block_idx][allele_idx]` → list of haplotype indices
 *   that correspond to this allele in the given block.
 *
 * Example:
 * ```
 * reverse_maps[b][a] = { h0, h1, h7 }
 * ```
 * means in block `b`, allele `a` corresponds to haplotypes 0, 1, and 7.
 *
 * @param[in] typed_only_reference_data Reference haplotype data partitioned
 *                                      into blocks, each containing allele
 *                                      cardinalities and a unique haplotype map.
 *
 * @return A nested vector structure containing reverse maps for each block.
 *
 * @complexity O(N), where N = total number of haplotype entries across all blocks.
 *
 * @note Each block in `typed_only_reference_data` must have consistent
 *       `cardinalities()` and `unique_map()` values.
 */
std::vector<std::vector<std::vector<std::size_t>>> generate_reverse_maps(const reduced_haplotypes& typed_only_reference_data);

/**
 * @brief Converts an old M3VCF file (v1/v2) to a newer VCF-like format (MVCFv3.0).
 *
 * This function reads an input M3VCF file, updates its headers and records,
 * and writes them to an output file in a modernized format compatible with
 * downstream tools. Optionally, a genetic map file may be provided to
 * annotate recombination positions.
 *
 * Conversion steps include:
 * - Parsing and updating VCF/M3VCF header lines.
 * - Ensuring presence of required `phasing` and `contig` headers.
 * - Adding INFO and FORMAT metadata fields required by MVCFv3.0.
 * - Reading sample IDs from the column header line.
 * - Deserializing haplotype blocks from the input file.
 * - Optionally annotating variants with centimorgan (cM) positions if a map file is provided.
 * - Serializing blocks to the output file in the new format.
 *
 * @param[in] input_path Path to the old M3VCF file (gzipped).
 * @param[in] output_path Path to the output file (can be `.bcf` or `.sav`).
 * @param[in] map_file_path Optional path to a genetic map file for cM annotation.
 *
 * @return True if the conversion completed successfully, false otherwise.
 *
 * @throws std::runtime_error If input or output files cannot be opened.
 *
 * @note
 * - Supports M3VCF v1 and v2 input.
 * - Output headers are modified to include MVCFv3.0 metadata.
 * - If `map_file_path` is non-empty, records will include cM positions.
 * - The function ensures `phasing` and `contig` headers exist.
 * - Sample IDs are extracted from the first non-header line of the input.
 *
 * @complexity
 * - Time: O(B × V), where B = number of blocks, V = number of variants per block.
 * - Memory: O(N), proportional to the size of haplotype data buffered during conversion.
 */
bool convert_old_m3vcf(const std::string& input_path, const std::string& output_path, const std::string& map_file_path = "");

/**
 * @brief Compress a haplotype reference panel into blocks and write to an output file.
 *
 * This function reads a phased reference panel from an input file (VCF/BCF/SAV format),
 * compresses haplotypes into blocks of unique haplotypes, and writes the result to
 * an output file in SAV/BCF format. It adaptively determines when to flush blocks
 * based on compression ratio and block size constraints.
 *
 * @param input_path       Path to the input reference panel file (VCF/BCF/SAV).
 * @param output_path      Path to the compressed output file (SAV/BCF).
 * @param min_block_size   Minimum number of variants required in a block before flushing.
 * @param max_block_size   Maximum number of variants allowed in a block before forcing flush.
 * @param slope_unit       Interval of variants used to check compression ratio slope.
 * @param map_file_path    Path to a genetic map file (currently unused, reserved for CM filling).
 *
 * @return true if compression and writing completed successfully, false otherwise.
 *
 * @note
 *  - Input file must contain fully phased genotypes (phasing header must not be "none" or "partial").
 *  - INFO and FORMAT fields required for compressed blocks are automatically added to headers.
 *  - Compression ratio is defined as:
 *    \f[
 *    CR = \frac{\text{expanded haplotype size} + 
 *    (\text{unique haplotype size} \times \text{variant size})}
 *    {\text{expanded haplotype size} \times \text{variant size}}
 *    \f]
 *
 * @see unique_haplotype_block, reference_site_info
 */
bool compress_reference_panel(const std::string& input_path, const std::string& output_path,
  std::size_t min_block_size = 10,
  std::size_t max_block_size = 0xFFFF, // max s1r block size minus 1 partition record
  std::size_t slope_unit = 10, 
  const std::string& map_file_path = "");


#endif // MINIMAC4_INPUT_PREP_HPP
