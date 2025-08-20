#ifndef MINIMAC4_HIDDEN_MARKOV_MODEL_HPP
#define MINIMAC4_HIDDEN_MARKOV_MODEL_HPP

#include "unique_haplotype.hpp"
#include "variant.hpp"

#include <array>
#include <cassert>
#include <numeric>

/**
 * @brief Stores full and leave-one-out (LOO) dosages for imputed variants.
 *
 * This class holds the imputed genotype probabilities for all variants and
 * haplotypes in the dataset. It provides storage for both standard dosages
 * and leave-one-out dosages, which are used for model validation and
 * cross-checking accuracy.
 */
class full_dosages_results
{
public:
  /**
   * @brief Matrix of imputed dosages.
   *
   * Each row corresponds to a variant, and each column corresponds to a
   * target haplotype or sample. Values are floating-point dosages in [0,1].
   */
  std::vector<std::vector<float>> dosages_;

  /**
   * @brief Matrix of leave-one-out (LOO) dosages.
   *
   * These dosages are computed excluding the current target haplotype from
   * the reference panel, useful for cross-validation of imputation.
   */
  std::vector<std::vector<float>> loo_dosages_;
public:
  /**
   * @brief Resizes the dosage matrices to the specified dimensions.
   *
   * This function resizes two internal 2D vectors: `dosages_` and `loo_dosages_`.
   * Each element is initialized to a sentinel value returned by
   * `savvy::typed_value::end_of_vector_value<float>()`.
   *
   * @param n_rows Number of rows for the main `dosages_` matrix.
   * @param n_loo_rows Number of rows for the leave-one-out `loo_dosages_` matrix.
   * @param n_columns Number of columns for both matrices.
   *
   * @details
   * After calling this function:
   * - `dosages_` will have dimensions `n_rows x n_columns`.
   * - `loo_dosages_` will have dimensions `n_loo_rows x n_columns`.
   * - All elements in both matrices are initialized to a sentinel value
   *   indicating the end of a vector, typically used for safe iteration or
   *   checking uninitialized entries.
   */
  void resize(std::size_t n_rows, std::size_t n_loo_rows, std::size_t n_columns)
  {
    dosages_.resize(n_rows, std::vector<float>(n_columns, savvy::typed_value::end_of_vector_value<float>()));
    loo_dosages_.resize(n_loo_rows, std::vector<float>(n_columns, savvy::typed_value::end_of_vector_value<float>()));
  }

  /**
   * @brief Fills all dosage matrices with the end-of-vector sentinel value.
   *
   * This function sets every element in the `dosages_` and `loo_dosages_`
   * 2D vectors to a sentinel value returned by
   * `savvy::typed_value::end_of_vector_value<float>()`.
   *
   * @details
   * This is useful for reinitializing matrices before recalculation or
   * ensuring that all entries are in a known "empty" state.
   * 
   * After calling this function:
   * - All elements in `dosages_` are set to the sentinel value.
   * - All elements in `loo_dosages_` are set to the sentinel value.
   */
  void fill_eov()
  {
    for (auto& v : dosages_)
      std::fill(v.begin(), v.end(),savvy::typed_value::end_of_vector_value<float>());

    for (auto& v : loo_dosages_)
      std::fill(v.begin(), v.end(),savvy::typed_value::end_of_vector_value<float>());
  }

  /**
   * @brief Returns the dimensions of the main dosages matrix.
   *
   * @return A `std::array` of size 2:
   *         - [0]: number of rows (`dosages_.size()`)
   *         - [1]: number of columns (`dosages_[0].size()` if not empty, 0 otherwise)
   */
    std::array<std::size_t, 2> dimensions() const { return {dosages_.size(), dosages_.empty() ? 0 : dosages_[0].size()}; }


  /**
   * @brief Returns the dimensions of the leave-one-out dosages matrix.
   *
   * @return A `std::array` of size 2:
   *         - [0]: number of rows (`loo_dosages_.size()`)
   *         - [1]: number of columns (`loo_dosages_[0].size()` if not empty, 0 otherwise)
   */
    std::array<std::size_t, 2> dimensions_loo() const { return {dosages_.size(), dosages_.empty() ? 0 : dosages_[0].size()}; }


  /**
   * @brief Accesses an element of the main dosages matrix (modifiable).
   *
   * @param i Row index.
   * @param j Column index.
   * @return Reference to the element at position (i, j) in `dosages_`.
   */
  float& dosage(std::size_t i, std::size_t j) { return dosages_[i][j]; }

  /**
   * @brief Accesses an element of the main dosages matrix (const version).
   *
   * @param i Row index.
   * @param j Column index.
   * @return Const reference to the element at position (i, j) in `dosages_`.
   */
  const float& dosage(std::size_t i, std::size_t j) const  { return dosages_[i][j]; }

  /**
   * @brief Accesses an element of the leave-one-out dosages matrix (modifiable).
   *
   * @param i Row index.
   * @param j Column index.
   * @return Reference to the element at position (i, j) in `loo_dosages_`.
   */
  float& loo_dosage(std::size_t i, std::size_t j) { return loo_dosages_[i][j]; }

  /**
   * @brief Accesses an element of the leave-one-out dosages matrix (const version).
   *
   * @param i Row index.
   * @param j Column index.
   * @return Const reference to the element at position (i, j) in `loo_dosages_`.
   */
  const float& loo_dosage(std::size_t i, std::size_t j) const { return loo_dosages_[i][j]; }
};

/**
 * @brief Implements a Hidden Markov Model for genotype imputation.
 *
 * This class performs multi-stage HMM-based imputation of genotype dosages
 * for target haplotypes using reference haplotype blocks. It maintains
 * forward and backward probability matrices, junction proportions, and
 * intermediate haplotype states for S1, S2, and S3 probability transformations.
 *
 * @details
 * - The model supports typed and untyped site imputation.
 * - Probabilities are normalized and adjusted for recombination, background error,
 *   and leave-one-out cross-validation.
 * - Implements forward traversal (`traverse_forward`) and backward traversal
 *   (`traverse_backward`) along haplotype blocks.
 */
class hidden_markov_model
{
private:
  /** Forward probability matrices per haplotype block. */
  std::deque<std::vector<std::vector<float>>> forward_probs_;

  /** Forward probabilities ignoring recombination. */
  std::deque<std::vector<std::vector<float>>> forward_norecom_probs_;

  /** Junction probability proportions per haplotype block. */
  std::vector<std::vector<float>> junction_prob_proportions_;

  /** Flags indicating if a precision jump occurred at a given site. */
  std::vector<bool> precision_jumps_;

  /** Minimum posterior probability threshold for accepting imputed states. */
  float prob_threshold_ = 0.01f;

  /** Threshold for S1 haplotype probability (used in multi-stage imputation). */
  float s1_prob_threshold_ = -1.f;

  /** Difference threshold for updating best haplotypes. */
  float diff_threshold_ = 0.01f;

  /** Background sequencing/genotyping error rate. */
  float background_error_ = 1e-5f;

  /** Recombination decay factor for distance-based probability adjustment. */
  double decay_ = 0.;

  /** Factor used to scale probabilities if they become too small. */
  static constexpr float jump_fix = 1e15f;

  /** Threshold for identifying underflow in probabilities. */
  static constexpr float jump_threshold = 1e-10f;

  /** Scaling factor for discretizing dosages. */
  const std::int16_t bin_scalar_ = 1000;

  /** Best S1 haplotypes (indices) for current variant site. */
  std::vector<std::uint32_t> best_s1_haps_;

  /** Best S2 haplotypes (indices) for current variant site. */
  std::vector<std::uint32_t> best_s2_haps_;

  /** Best S3 haplotypes (indices) for current variant site. */
  std::vector<std::uint32_t> best_s3_haps_;

  /** Posterior probabilities corresponding to best S1 haplotypes. */
  std::vector<float> best_s1_probs_;

  /** Posterior probabilities corresponding to best S2 haplotypes. */
  std::vector<float> best_s2_probs_;

  /** S2 posterior probabilities used in intermediate computations. */
  std::vector<float> s2_probs_;

  /** Cardinalities of S2 haplotypes (number of expanded haplotypes per unique haplotype). */
  std::vector<std::size_t> s2_cardinalities_;

  /** Posterior probabilities corresponding to best S3 haplotypes. */
  std::vector<float> best_s3_probs_;

public:
  /**
   * @brief Constructs a Hidden Markov Model with specified parameters.
   *
   * @param s3_prob_threshold Threshold probability for S3 state.
   * @param s1_prob_threshold Threshold probability for S1 state.
   * @param diff_threshold Minimum difference required between probabilities to
   *                       make a confident state call.
   * @param background_error Expected background error rate.
   * @param decay Decay factor controlling the influence of previous states.
   *
   * @details
   * This constructor initializes the internal HMM parameters. These thresholds
   * and the decay factor influence the model's sensitivity to differences in
   * observed probabilities and determine how the hidden states are inferred.
   */
  hidden_markov_model(float s3_prob_threshold, float s1_prob_threshold, float diff_threshold, float background_error, float decay);

  /**
   * @brief Performs a forward traversal over reference haplotypes for a given target haplotype.
   *
   * This function implements the forward algorithm of a Hidden Markov Model (HMM)
   * to compute the probability of observing the target variants given the reference
   * haplotypes. It accounts for recombination events and optional leave-one-out 
   * calculations.
   *
   * @param ref_haps A deque of reference haplotype blocks (`unique_haplotype_block`).
   *                 Each block contains multiple haplotypes and variant positions.
   * @param tar_variants A vector of target variants (`target_variant`) for which
   *                    the forward probabilities are computed.
   * @param hap_idx The index of the target haplotype to traverse within `tar_variants`.
   *
   * @details
   * The function performs the following steps:
   * 1. Initializes internal forward probability matrices (`forward_probs_`, `forward_norecom_probs_`)
   *    and junction probabilities (`junction_prob_proportions_`).
   * 2. For the first reference block, likelihoods are initialized using `initialize_likelihoods`.
   * 3. For subsequent blocks:
   *    - Computes junction probabilities based on previous block probabilities and
   *      recombination proportions.
   *    - Normalizes junction probabilities to ensure they sum to 1.
   *    - Updates `precision_jumps_` by transposing probabilities to the next position.
   * 4. Iterates over all variants in each block:
   *    - Updates probabilities conditioned on observed target genotypes using `condition`.
   *    - Computes forward probabilities for the next position using `transpose`.
   *
   * @note
   * - Assertions ensure probability values remain within valid ranges [0, 1].
   * - The function supports handling missing genotypes (`observed < 0`) by skipping
   *   conditioning for those positions.
   */
  void traverse_forward(const std::deque<unique_haplotype_block>& ref_haps,
    const std::vector<target_variant>& tar_variant,
    std::size_t hap_idx);

  /**
   * @brief Performs a backward traversal over reference haplotypes to compute posterior probabilities.
   *
   * This function implements the backward algorithm of a Hidden Markov Model (HMM)
   * for a given target haplotype. It computes probabilities by combining forward
   * probabilities (from `traverse_forward`) with backward likelihoods, handling
   * recombination events, junction proportions, and leave-one-out considerations.
   *
   * @param ref_haps A deque of reference haplotype blocks (`unique_haplotype_block`).
   *                 Each block contains multiple haplotypes and variant positions.
   * @param tar_variants A vector of target variants (`target_variant`) to traverse backward.
   * @param hap_idx Index of the target haplotype to traverse within `tar_variants`.
   * @param out_idx Output index used for storing results in `full_dosages_results`.
   * @param reverse_maps Maps from expanded to unique haplotypes for each block.
   * @param output Reference to a `full_dosages_results` structure where computed
   *               posterior dosages are stored.
   * @param full_reference_data Reduced haplotype reference data needed for imputation.
   *
   * @details
   * The function performs the following steps:
   * 1. Initializes backward probability vectors (`backward`, `backward_norecom`) and
   *    junction proportions for the last block.
   * 2. Iterates backward over all reference blocks:
   *    - Computes junction probabilities using the backward pass.
   *    - Normalizes probabilities to ensure valid distributions.
   *    - Updates constants used for imputation with combined forward and backward probabilities.
   *    - Calls `impute` to compute posterior dosages for each variant position.
   *    - Conditions backward probabilities on observed target genotypes if available.
   * 3. Handles recombination corrections and precision jumps using `transpose`.
   *
   * @note
   * - Assertions ensure probability values remain within valid ranges and that
   *   the global index reaches -1 at the end.
   * - Observed genotypes with negative values are treated as missing and skipped
   *   during conditioning.
   * - This method relies on data structures populated during `traverse_forward`.
   */
  void traverse_backward(const std::deque<unique_haplotype_block>& ref_haps,
    const std::vector<target_variant>& tar_variant,
    std::size_t hap_idx,
    std::size_t out_idx,
    const std::vector<std::vector<std::vector<std::size_t>>>& reverse_maps,
    full_dosages_results& output,
    const reduced_haplotypes& full_reference_data);
private:
  /**
   * @brief Updates forward or backward probabilities conditioned on an observed genotype.
   *
   * This function adjusts the probabilities of each haplotype given the observed
   * genotype at a specific variant, accounting for genotyping error and background error.
   *
   * @param probs Vector of probabilities to be updated (e.g., forward or backward probabilities).
   * @param probs_norecom Vector of probabilities without recombination to be updated.
   * @param template_haps Vector of template haplotype genotypes at the current variant.
   * @param observed Observed genotype of the target haplotype at this variant.
   *                 A non-negative integer represents a valid genotype.
   * @param err Genotyping error probability.
   * @param af Allele frequency of the observed variant.
   *
   * @details
   * - Computes a probability of random mismatch (`prandom`) based on error and allele frequency.
   * - Computes the match probability (`pmatch`) for haplotypes matching the observed genotype.
   * - Multiplies each probability by either `pmatch` or `prandom` depending on whether
   *   the template haplotype matches the observation.
   * - Ensures that all resulting probabilities are non-negative.
   *
   * @note
   * - This function is typically called during the forward or backward traversal
   *   to condition probabilities on observed data.
   * - Assertions ensure probabilities remain non-negative after conditioning.
   */
  void condition(std::vector<float>& probs, std::vector<float>& probs_norecom, const std::vector<std::int8_t>& template_haps, std::int8_t observed, float err, float freq);
  
  /**
   * @brief Transposes probability vectors to the next variant position, accounting for recombination.
   *
   * This function updates the forward or backward probability vectors for the next position
   * in the Hidden Markov Model, distributing probabilities according to recombination rates
   * and haplotype cardinalities.
   *
   * @param from Input probability vector at the current position.
   * @param to Output probability vector for the next position.
   * @param from_norecom Input probability vector ignoring recombination.
   * @param to_norecom Output probability vector for the next position ignoring recombination.
   * @param uniq_cardinalities Vector of unique haplotype counts for each template.
   * @param recom Recombination probability between the current and next positions.
   * @param n_templates Number of templates (expanded haplotypes) at this position.
   * 
   * @return `true` if a probability "jump" correction was applied to prevent underflow; otherwise `false`.
   *
   * @details
   * - Computes the sum of the current probabilities.
   * - Applies a scaling factor if probabilities are too small (less than `jump_threshold`),
   *   using `jump_fix` to avoid numerical underflow.
   * - Updates `to` by combining the recombination-adjusted fraction and the scaled sum
   *   weighted by haplotype cardinalities.
   * - Updates `to_norecom` by scaling the non-recombination probabilities.
   * - Ensures all probabilities remain within a safe range and non-negative.
   *
   * @note
   * - Assertions check for probability validity after computation.
   * - This method is used internally in forward/backward traversals for HMM computations.
   */
  bool transpose(const std::vector<float>& from, std::vector<float>& to, const std::vector<float>& from_norecom, std::vector<float>& to_norecom, const std::vector<std::size_t>& uniq_cardinalities, double recom, std::size_t n_templates);

  /**
   * @brief Computes posterior dosages and probabilities for a single typed variant site.
   *
   * This function calculates the posterior probability that a target haplotype
   * carries the alternative allele at a given variant site, based on forward
   * and backward HMM probabilities, recombination-adjusted junction proportions,
   * and observed genotypes.
   *
   * @param prob_sum Reference to the sum of probabilities for normalization.
   * @param prev_best_hap Reference to the index of the previously best haplotype.
   *                      May be updated if a single haplotype exceeds the probability threshold.
   * @param left_probs Forward probabilities at the current position.
   * @param right_probs Backward probabilities at the current position.
   * @param left_probs_norecom Forward probabilities ignoring recombination.
   * @param right_probs_norecom Backward probabilities ignoring recombination.
   * @param left_junction_proportions Junction proportions from forward pass.
   * @param right_junction_proportions Junction proportions from backward pass.
   * @param constants Precomputed constants combining forward and backward contributions.
   * @param reverse_map Maps from expanded haplotypes to unique haplotypes.
   * @param template_haps Template haplotypes (0 or 1 for reference/alternate allele).
   * @param observed Observed genotype for the target haplotype.
   * @param err Genotyping error probability.
   * @param af Allele frequency of the variant.
   * @param best_unique_haps Output vector of indices for haplotypes exceeding the probability threshold.
   * @param best_unique_probs Output vector of probabilities for the corresponding haplotypes.
   * @param dose Output posterior dosage (0-1) for the typed site.
   * @param loo_dose Output leave-one-out dosage for the typed site.
   *
   * @details
   * - Checks if the previously best haplotype exceeds the threshold and assigns the dosage directly.
   * - Otherwise, computes posterior probabilities for all haplotypes, summing contributions for reference
   *   and alternative alleles.
   * - Normalizes the dosage and optionally bins it for numerical stability.
   * - Computes leave-one-out dosages adjusting for observed genotypes and genotyping error.
   * - Updates `prev_best_hap` to optimize repeated computations when a single haplotype dominates the probability.
   *
   * @note
   * - Probabilities are clipped and binned to avoid floating-point precision issues.
   * - Assertions ensure that intermediate probabilities are valid and positive.
   * - This function is a core step in computing posterior dosages for HMM-based genotype imputation.
   */
  void impute_typed_site(double& prob_sum, std::size_t& prev_best_hap,
    const std::vector<float>& left_probs,
    const std::vector<float>& right_probs,
    const std::vector<float>& left_probs_norecom,
    const std::vector<float>& right_probs_norecom,
    const std::vector<float>& left_junction_proportions,
    const std::vector<float>& right_junction_proportions,
    const std::vector<float>& constants,
    const std::vector<std::vector<std::size_t>>& reverse_map,
    const std::vector<std::int8_t>& template_haps,
    std::int8_t observed, float err, float af,
    std::vector<std::uint32_t>& best_uniq_haps, std::vector<float>& best_uniq_probs, float& dose, float& loo_dose);
  
  /**
   * @brief Imputes dosages for a target variant site using HMM forward/backward probabilities.
   *
   * This function computes posterior dosages and leave-one-out dosages for a single
   * variant site by combining forward and backward probabilities with junction proportions
   * and template haplotypes. It also updates intermediate HMM states for S3 → S1 → S2
   * probability transformations used in multi-stage imputation.
   *
   * @param prob_sum Reference to the sum of probabilities for normalization.
   * @param prev_best_typed_hap Reference to the previously best haplotype index for typed sites.
   *                            May be updated for optimization.
   * @param left_probs Forward probabilities at the current variant site.
   * @param right_probs Backward probabilities at the current variant site.
   * @param left_probs_norecom Forward probabilities ignoring recombination.
   * @param right_probs_norecom Backward probabilities ignoring recombination.
   * @param left_junction_proportions Junction proportions from forward pass.
   * @param right_junction_proportions Junction proportions from backward pass.
   * @param constants Precomputed constants combining forward and backward contributions.
   * @param uniq_map Map from expanded haplotypes to unique haplotypes.
   * @param reverse_map Reverse mapping from unique to expanded haplotypes.
   * @param template_haps Template haplotypes (0/1) at the current site.
   * @param tar_variants Vector of target variants.
   * @param row Index of the current variant in `tar_variants`.
   * @param column Index of the haplotype column for the target genotype.
   * @param out_column Column index in the output dosage matrix.
   * @param output Reference to `full_dosages_results` for storing imputed dosages.
   * @param full_ref_ritr Iterator to the current position in reduced reference haplotypes.
   * @param full_ref_rend Iterator to the end of the reduced reference haplotypes.
   * @param prev_block_idx Reference to previous block index for probability transformations.
   *
   * @details
   * - Calls `impute_typed_site` to compute typed site posterior probabilities and dosages.
   * - Compares results to previously stored S3 haplotypes and updates if differences exceed
   *   the threshold.
   * - Performs S3 → S1 → S2 probability transformations for multi-level imputation.
   * - Iterates over full reference haplotypes and updates output dosages:
   *   - For typed sites, assigns direct dosages.
   *   - For untyped sites, computes posterior probability of alternative alleles,
   *     adjusting for cardinalities, missing data, and recombination decay.
   * - Uses `bin_scalar_` to discretize probabilities for numerical stability.
   *
   * @note
   * - The method ensures probabilities remain within [0,1] and uses assertions to validate intermediate values.
   * - Handles first/last variant decay adjustments using recombination estimates.
   * - This method is a key step in backward imputation traversal, bridging typed and untyped sites.
   */
  void impute(double& prob_sum, std::size_t& prev_best_expanded_hap,
    const std::vector<float>& left_probs,
    const std::vector<float>& right_probs,
    const std::vector<float>& left_probs_norecom,
    const std::vector<float>& right_probs_norecom,
    const std::vector<float>& left_junction_proportions,
    const std::vector<float>& right_junction_proportions,
    const std::vector<float>& constants,
    const std::vector<std::int64_t>& uniq_map,
    const std::vector<std::vector<std::size_t>>& reverse_map,
    const std::vector<std::int8_t>& template_haps,
    const std::vector<target_variant>& tar_variants,
    std::size_t row, std::size_t column, std::size_t out_column,
    full_dosages_results& output,
    reduced_haplotypes::iterator& full_ref_ritr, const reduced_haplotypes::iterator& full_ref_rend,
    std::size_t& prev_block_idx);

  void initialize_likelihoods(std::vector<float>& probs, std::vector<float>& probs_norecom, std::vector<float>& proportions, const unique_haplotype_block& ref_block);

  void s3_to_s1_probs(
    const std::vector<float>& left_probs, const std::vector<float>& right_probs,
    const std::vector<float>& left_probs_norecom, const std::vector<float>& right_probs_norecom,
    const std::vector<float>& left_junction_proportions, const std::vector<float>& right_junction_proportions,
    const std::vector<std::vector<std::size_t>>& s3_reverse_map, double prob_sum);
  void s1_to_s2_probs(std::vector<std::size_t>& cardinalities, const std::vector<std::int64_t>& uniq_map, std::size_t s2_size);
};

#endif // MINIMAC4_HIDDEN_MARKOV_MODEL_HPP