#ifndef MINIMAC4_UNIQUE_HAPLOTYPE_HPP
#define MINIMAC4_UNIQUE_HAPLOTYPE_HPP

#include "variant.hpp"
#include "recombination.hpp"

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <limits>
#include <cassert>

/**
 * @class unique_haplotype_block
 * @brief Represents a block of unique haplotypes and their variants.
 *
 * This class stores haplotype information in a compressed form, mapping individual haplotypes
 * to unique columns of alleles, tracking allele counts (cardinalities), and storing variant details.
 */
class unique_haplotype_block
{
private:
  /** 
   * @brief Maps each haplotype to a unique allele column index.
   *
   * Each element corresponds to a haplotype and stores the index of the
   * unique allele column it belongs to. End-of-vector sentinel values
   * indicate missing or invalid haplotypes.
   */
  std::vector<std::int64_t> unique_map_;

  /** 
   * @brief Cardinalities of unique haplotype columns.
   *
   * Stores the number of haplotypes assigned to each unique allele column.
   * The size of this vector matches the number of unique haplotypes.
   */
  std::vector<std::size_t> cardinalities_;

  /** 
   * @brief List of variants for this haplotype block.
   *
   * Each element contains detailed information about a variant, including
   * position, reference and alternate alleles, recombination rates, and
   * haplotype genotypes.
   */
  std::vector<reference_variant> variants_;
public:
  /**
   * @brief Compress and map haplotype alleles for a new variant into the block.
   *
   * This function updates the block of unique haplotypes by incorporating a new
   * variant, compressing redundant haplotypes, and maintaining a mapping from
   * input alleles to unique haplotype indices.
   *
   * The first call initializes the block with the given variant. Subsequent
   * calls attempt to map the alleles to existing haplotype columns or create
   * new ones if necessary.
   *
   * @param site_info  Reference site metadata (chromosome, position, ID, alleles, error, recombination rate, cM).
   * @param alleles    Vector of observed alleles (per haplotype).
   *
   * @return true if compression and mapping succeeded, false otherwise.
   *
   * @details
   * - If this is the **first variant**, a unique mapping is initialized:
   *   - Each unique allele is assigned a haplotype column.
   *   - A `unique_map_` is built from haplotype indices to allele indices.
   *   - `cardinalities_` tracks the number of haplotypes assigned per column.
   *
   * - For **subsequent variants**:
   *   - Alleles must match the size of `unique_map_`, otherwise the function fails.
   *   - Each haplotype is mapped back to its column via `unique_map_`.
   *   - If the allele matches the expected value, it is stored directly.
   *   - If it mismatches:
   *     - Check if it matches a newly created column.
   *     - If not, a **new haplotype column** is created and propagated across all variants.
   *
   * - Special handling:
   *   - `savvy::typed_value::is_end_of_vector()` is used to mark missing/invalid alleles.
   *   - Ploidy mismatches are detected and reported as errors.
   *   - `cardinalities_` is updated to reflect haplotype counts per column.
   *
   * @note
   * - Ensures that the total haplotype count matches across all variants.
   * - Uses `assert()` checks to verify internal consistency.
   * - Returns `false` if:
   *   - The alleles vector is empty.
   *   - A sample ploidy mismatch is detected.
   *   - The haplotype size differs from the expected mapping.
   *
   * Example use:
   * @code
   * unique_haplotype_block block;
   * reference_site_info site{...};
   * std::vector<std::int8_t> alleles = {0, 1, 0, 1};
   * bool ok = block.compress_variant(site, alleles);
   * @endcode
   */
  bool compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles);

  /**
   * @brief Get the list of compressed reference variants in this block.
   * @return Const reference to the vector of reference variants.
   */
  const std::vector<reference_variant>& variants() const { return variants_; }

  /**
   * @brief Get the mapping from haplotypes to unique columns.
   *
   * Each entry corresponds to an original haplotype index and indicates
   * which compressed column it maps to.
   *
   * @return Const reference to the unique haplotype mapping vector.
   */
  const std::vector<std::int64_t>& unique_map() const { return unique_map_; }

  /**
   * @brief Get the number of original haplotypes after expansion.
   *
   * This corresponds to the length of the `unique_map_` vector.
   *
   * @return Expanded haplotype count.
   */
  std::size_t expanded_haplotype_size() const { return unique_map_.size(); }

  /**
   * @brief Get the number of unique haplotypes represented in the block.
   *
   * This corresponds to the number of columns in the genotype matrix.
   * If the block is empty, returns 0.
   *
   * @return Unique haplotype count.
   */
  std::size_t unique_haplotype_size() const { return variants_.empty() ? 0 : variants_[0].gt.size(); }

  /**
   * @brief Get the number of variants compressed in this block.
   * @return Variant count.
   */
  std::size_t variant_size() const { return variants_.size(); }

  /**
   * @brief Get the haplotype cardinalities across unique columns.
   *
   * Each entry represents how many haplotypes map to the corresponding
   * unique haplotype column.
   *
   * @return Const reference to the cardinalities vector.
   */
  const std::vector<std::size_t>& cardinalities() const { return cardinalities_; }

  /**
   * @brief Clears all data stored in the haplotype block.
   *
   * This method resets the internal state of the haplotype block by removing all
   * stored variants, haplotype mappings, and cardinality counts.
   *
   * Specifically:
   * - `variants_` is cleared, removing all variant records.
   * - `unique_map_` is cleared, removing all haplotype mapping indices.
   * - `cardinalities_` is cleared, resetting the frequency counts of haplotypes.
   *
   * After calling this method, the haplotype block is in an empty state
   * and must be repopulated (e.g., via `deserialize()`).
   */
  void clear();

  /**
   * @brief Trims variants outside a specified genomic range.
   *
   * This function removes variants from the haplotype block that fall
   * outside the interval `[min_pos, max_pos]`.
   *
   * Behavior:
   * - If all variants are outside the range, the block is cleared entirely.
   * - Otherwise, variants before `min_pos` or after `max_pos` are erased,
   *   while the remaining variants are preserved.
   *
   * @param min_pos The minimum genomic position (inclusive).
   * @param max_pos The maximum genomic position (inclusive).
   *
   * @note This function only modifies the `variants_` container. The mappings
   *       (`unique_map_`) and cardinalities remain unchanged, so care must
   *       be taken when trimming after compression.
   */
  void trim(std::size_t min_pos, std::size_t max_pos);

  /**
   * @brief Removes the last variant from the haplotype block.
   *
   * This function erases the most recently added variant in the `variants_` vector.
   *
   * @note This only modifies the `variants_` container. It does not update
   *       `unique_map_` or `cardinalities_`, so use with caution if the block
   *       has been compressed or mapped.
   */
  void pop_variant();

  /**
   * @brief Fills the centimorgan (cM) values for all variants in the haplotype block.
   *
   * This function iterates over all variants in the `variants_` vector and
   * sets each variant's `cm` field by interpolating its genetic position
   * using the provided `genetic_map_file`.
   *
   * @param map_file Reference to a `genetic_map_file` object used for interpolation.
   */
  void fill_cm(genetic_map_file& map_file);

  /**
   * @brief Fills missing centimorgan (cM) values for variants using recombination probabilities.
   *
   * This function iterates over all variants in the `variants_` vector. For any variant
   * with a NaN `cm` value, it sets `cm` to the provided `start_cm`. If the variant has a
   * valid `recom` value, `start_cm` is incremented using the recombination probability
   * converted to cM via `recombination::switch_prob_to_cm()`.
   *
   * @param start_cm Reference to a double representing the starting centimorgan value.
   *                 This value will be updated as variants are processed.
   */
  void fill_cm_from_recom(double& start_cm);

  /**
   * @brief Deserialize a unique haplotype block from an input stream (m3vcf format).
   *
   * This function reads and parses a block of data from a multiple-phased VCF (m3vcf) file
   * into the `unique_haplotype_block` object. It populates internal structures such as
   * the unique haplotype map, variant information, and cardinalities.
   *
   * @param[in,out] is            Input stream containing m3vcf block data. Stream state
   *                              is updated depending on success or error.
   * @param[in]     m3vcf_version Format version of the m3vcf file (1 or 2).
   * @param[in]     n_haplotypes  Total number of haplotypes expected in this block.
   *
   * @return `true` if deserialization was successful, `false` otherwise.
   *
   * @note On failure, the method clears all internal data structures, marks the stream
   *       with `std::ios::badbit`, and writes an error message to `std::cerr`.
   *
   * ## Parsing Steps
   * 1. Read the header line of the block, extract metadata:
   *    - `VARIANTS=<N>` → number of variants in the block.
   *    - `REPS=<M>` → number of unique haplotype representatives.
   * 2. Parse haplotype indices (`unique_map_`) from genotype columns:
   *    - For version 1: one haplotype index per column.
   *    - For version 2: paired haplotype indices separated by '|'.
   * 3. Validate that the number of parsed haplotypes matches `n_haplotypes`.
   * 4. Build the `cardinalities_` vector (frequency of each unique haplotype).
   * 5. Read `n_variants` subsequent lines, parsing:
   *    - Chromosome, position, ID, reference/alternate alleles.
   *    - INFO fields (extracts `ERR` and `RECOM` values).
   *    - Genotypes:
   *      - For version 2: encoded as run-length style offsets (`cols[8]`).
   *      - For version 1: raw 0/1 vector of length `n_reps`.
   *    - Allele counts (`ac`) computed using `cardinalities_`.
   *
   * ## Error Conditions
   * - Invalid format (wrong number of columns, invalid integer conversions, etc.).
   * - Number of haplotypes read does not match `n_haplotypes`.
   * - Genotype column inconsistent with expected `n_reps`.
   * - Truncated input (fewer than expected variant lines).
   *
   * In any of these cases, the function:
   * - Clears internal state (`clear()`).
   * - Marks the stream as bad (`is.setstate(is.rdstate() | std::ios::badbit)`).
   * - Prints an error message to `std::cerr`.
   */
  bool deserialize(std::istream& is, int m3vcf_version, std::size_t n_haplotypes);

  /**
   * @brief Deserializes a unique haplotype block from a SAVVY input file and variant.
   *
   * This function reads haplotype block data from the given `input_file` using
   * the `var` object. It populates the `variants_`, `unique_map_`, and `cardinalities_`
   * of the `unique_haplotype_block`. The function stops reading when a "<BLOCK>" alt
   * allele is encountered or the input ends.
   *
   * @param input_file Reference to a `savvy::reader` object representing the input file.
   * @param var Reference to a `savvy::variant` object used to read variant data.
   * @return int Returns:
   *   - `-1` on I/O errors or invalid data,
   *   - `0` if input is empty or not ready,
   *   - `n+1` where `n` is the number of variants successfully deserialized.
   *
   * @note The `gt` vector of each variant is expected to match the size of `cardinalities_`.
   *       The allele count `ac` is computed using the `cardinalities_`.
   *       The first "<BLOCK>" allele encountered marks the end of this haplotype block.
   */
  int deserialize(savvy::reader& input_file, savvy::variant& var); // return <0 is error, 0 is EOF, >0 good

  /**
   * @brief Serializes the unique haplotype block to a SAVVY output file.
   *
   * This function writes the current haplotype block, including its variants,
   * unique haplotype mapping, and cardinalities, to the specified `output_file`.
   * The first variant in the block is treated as a marker with "<BLOCK>" alt allele.
   * Each variant stores allele counts (AC/AN), error probabilities, recombination rates,
   * centimorgan positions, and haplotype genotypes (UHA).
   *
   * @param output_file Reference to a `savvy::writer` object representing the output file.
   * @return bool Returns `true` if serialization succeeded, `false` if the block is empty
   *              or on I/O errors.
   *
   * @note The function sets the block size in `output_file` to align zstd compression blocks
   *       with m3vcf haplotype blocks.
   */
  bool serialize(savvy::writer& output_file);

  /**
 * @brief Removes "end-of-vector" markers from the unique haplotype map.
 *
 * This function scans the `unique_map_` vector and removes any elements
 * that represent the end-of-vector (EOV) sentinel value used to indicate
 * missing or invalid haplotypes.
 *
 * After this operation, the `unique_map_` will contain only valid haplotype indices,
 * and its size is reduced by the number of EOV entries removed.
 */
  void remove_eov();
};


/**
 * @class reduced_haplotypes
 * @brief Represents a collection of haplotype blocks with reduced storage.
 *
 * This class manages a sequence of unique_haplotype_block objects, allowing
 * efficient storage and access to haplotypes across multiple genomic variants.
 * It maintains block offsets, variant counts, and supports compression of variants.
 *
 * @note Iterators traverse variants across all blocks in order.
 */
class reduced_haplotypes
{
private:
  std::vector<std::size_t> block_offsets_;
  std::deque<unique_haplotype_block> blocks_;
  std::size_t variant_count_ = 0;
  std::size_t min_block_size_ = 1;
  std::size_t max_block_size_ = std::numeric_limits<std::size_t>::max();
  bool flush_block_ = true;
public:

  /**
   * @class reduced_haplotypes::iterator
   * @brief Iterator for traversing variants within reduced_haplotypes.
   *
   * Provides forward and backward traversal across haplotype blocks.
   * Supports standard iterator operations (*, ->, ++, --) and index access.
   *
   * @note The iterator handles crossing block boundaries automatically.
   */
  class iterator
  {
  private:
    const reduced_haplotypes* parent_;
    std::size_t block_idx_;
    std::size_t variant_idx_;
  public:
    iterator(const reduced_haplotypes& parent, std::size_t block_idx, std::size_t variant_idx) :
      parent_(&parent),
      block_idx_(block_idx),
      variant_idx_(variant_idx)
    {
    }

    iterator& operator++()
    {
      ++variant_idx_;
      assert(block_idx_ < parent_->blocks_.size());
      if (variant_idx_ >= parent_->blocks_[block_idx_].variant_size())
      {
        ++block_idx_;
        variant_idx_ = 0;
      }
      return *this;
    }

    iterator operator++(int)
    {
      iterator ret(*this);
      ++(*this);
      return ret;
    }

    iterator& operator--()
    {
      if (variant_idx_ == 0)
      {
        --block_idx_;
        if (block_idx_ < parent_->blocks_.size())
        {
          (*this) = iterator(*parent_, block_idx_, parent_->blocks_[block_idx_].variant_size());
          --(*this);
        }
      }
      else
      {
        --variant_idx_;
      }

      return *this;
    }

    iterator operator--(int)
    {
      iterator ret(*this);
      --(*this);
      return ret;
    }

    const reference_variant& operator*() const
    {
      return parent_->blocks_[block_idx_].variants()[variant_idx_];
    }

    const reference_variant* operator->() const
    {
      return &(parent_->blocks_[block_idx_].variants()[variant_idx_]);
    }

    std::size_t block_idx() const { return block_idx_; }
    std::size_t block_local_idx() const { return variant_idx_; }
    std::size_t global_idx() const { return parent_->block_offsets_[block_idx_] + variant_idx_; }
    const std::vector<std::int64_t>& unique_map() const { return parent_->blocks_[block_idx_].unique_map(); }
    const std::vector<std::size_t>& cardinalities() const { return parent_->blocks_[block_idx_].cardinalities(); }
  };

  iterator begin() const { return iterator(*this, 0, 0); }
  iterator end() const { return iterator(*this, blocks_.size(), 0); }

  reduced_haplotypes() {}

  /**
   * @brief Constructs a reduced_haplotypes object with specified block size limits.
   *
   * This constructor initializes the minimum and maximum haplotype block sizes.
   * It ensures that both sizes are at least 1.
   *
   * @param min_block_size Minimum allowed size for a haplotype block.
   * @param max_block_size Maximum allowed size for a haplotype block.
   */
  reduced_haplotypes(std::size_t min_block_size, std::size_t max_block_size);

  /**
   * @brief Compresses a variant into the current reduced haplotype block.
   *
   * This function attempts to add a variant to the last haplotype block. If the 
   * compression ratio improves or the flush flag is set, a new block is started.
   *
   * @param site_info Reference information for the variant (chromosome, position, alleles, etc.).
   * @param alleles Vector of alleles corresponding to each haplotype.
   * @param flush_block If true, forces the current block to flush and start a new block.
   * @return True if the variant was successfully compressed into the block, false otherwise.
   */
  bool compress_variant(const reference_site_info& site_info, const std::vector<std::int8_t>& alleles, bool flush_block = false);
  
  /**
   * @brief Appends a new haplotype block to the collection of reduced haplotype blocks.
   *
   * If the first variant of the new block is identical to the last variant of the
   * previous block (same position, reference, and alternate alleles), the duplicate
   * variant from the previous block is removed to avoid redundancy.
   *
   * @param block The unique haplotype block to append.
   */
  void append_block(const unique_haplotype_block& block);

  /**
   * @brief Fills the centimorgan (cM) values for all variants in all blocks
   *        using the provided genetic map.
   *
   * @param map_file The genetic map file used to interpolate centimorgan values.
   */
  void fill_cm(genetic_map_file& map_file);

  /**
   * @brief Calculates the overall compression ratio of all haplotype blocks.
   *
   * The compression ratio is defined as the ratio between the sum of expanded haplotype sizes
   * plus the product of unique haplotype size and variant count, divided by the total expanded
   * haplotype size multiplied by the number of variants across all blocks.
   *
   * @return The compression ratio as a float.
   */
  float compression_ratio() const;

  /**
   * @brief Accesses the stored haplotype blocks.
   *
   * @return A constant reference to the deque of unique_haplotype_block objects.
   */
  const std::deque<unique_haplotype_block>& blocks() const { return blocks_; }

  /**
   * @brief Returns the total number of variants across all blocks.
   *
   * @return Total variant count as std::size_t.
   */
  std::size_t variant_size() const { return variant_count_; }
};

/**
 * @brief Compares two iterators for equality.
 *
 * Two iterators are considered equal if they refer to the same block index
 * and the same local index within that block.
 *
 * @param lhs Left-hand side iterator.
 * @param rhs Right-hand side iterator.
 * @return True if iterators are equal, false otherwise.
 */
inline bool operator==(const reduced_haplotypes::iterator& lhs, const reduced_haplotypes::iterator& rhs)
{
  return lhs.block_idx() == rhs.block_idx() && lhs.block_local_idx() == rhs.block_local_idx();
}

/**
 * @brief Compares two iterators for inequality.
 *
 * This is the negation of operator==.
 *
 * @param lhs Left-hand side iterator.
 * @param rhs Right-hand side iterator.
 * @return True if iterators are not equal, false otherwise.
 */
inline bool operator!=(const reduced_haplotypes::iterator& lhs, const reduced_haplotypes::iterator& rhs)
{
  return !(lhs == rhs);
}



#endif // MINIMAC4_UNIQUE_HAPLOTYPE_HPP