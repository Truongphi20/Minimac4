#ifndef MINIMAC4_RECOMBINATION_HPP
#define MINIMAC4_RECOMBINATION_HPP

#include "variant.hpp"

#include <shrinkwrap/istream.hpp>

#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>

/**
 * @class recombination
 * @brief Provides utilities for handling recombination rates and genetic maps.
 *
 * The `recombination` class is responsible for parsing genetic map files,
 * interpolating recombination rates across target variant sites, and converting
 * between genetic map units (centiMorgans) and switch probabilities used in
 * haplotype inference.
 *
 * ### Responsibilities
 * - Parse and align genetic map files to variant sites.
 * - Interpolate recombination map values for variants without direct entries.
 * - Convert between genetic distances (cM) and recombination/switch probabilities
 *   using models such as Haldane’s map function.
 *
 * ### Key Components
 * - `struct map_file_line`: Represents a single line in the map file,
 *   including chromosome, base-pair position, and map value.
 * - `parse_map_file()`: Aligns recombination map data with target variants
 *   and computes switch probabilities.
 * - `read_entry()`: Reads a line from the map file (supports both old and new formats).
 * - Utility functions:
 *   - `haldane()` and `haldane_inverse()` for recombination ↔ cM conversions.
 *   - `cm_to_switch_prob()` and `switch_prob_to_cm()` for cM ↔ switch probability conversions.
 *
 * ### Notes
 * - The last recombination probability in the sequence is always set to zero
 *   to ensure proper backward traversal during haplotype inference.
 * - Supports both the old PLINK-style map format and new three-column format.
 */
class recombination
{
public:
  //static constexpr float recom_min = 1e-5f;
//private:

  /**
   * @brief Representation of a single line from a genetic map file.
   *
   * Each line in the map file provides information about a genomic position
   * and its corresponding genetic map value (in centiMorgans).
   */
  struct map_file_line
  {
   /** @brief Chromosome identifier (e.g., "1", "chrX"). */
  std::string chrom;

  /** @brief Base-pair position on the chromosome. */
  std::size_t pos = 0;

  /** @brief Genetic map value at the position (in centiMorgans). */
  float map_value = 0.;
  };
public:
  /**
   * @brief Parse and align a genetic recombination map file to a set of target variants.
   *
   * This function reads a genetic map file and assigns recombination rates (in cM) 
   * to each target variant in the provided vector. It aligns the genetic map to 
   * the given chromosome, interpolates recombination values between map positions, 
   * and converts the aligned genetic map into switch probabilities for HMM traversal.
   *
   * The recombination probability at each site is computed based on the genetic 
   * distance between consecutive variants and is adjusted to ensure non-negative 
   * values, respecting the provided minimum recombination threshold.
   *
   * @param map_file_path Path to the genetic map file. Must contain either 
   *        standard format or a header-prefixed "new" format (three tab-delimited columns).
   * @param sites Vector of target variants (with chromosome, position, and recombination fields).
   *        The recombination field (`recom`) is updated in-place.
   * @param recom_min Minimum recombination probability to enforce for each site.
   *
   * @return true if parsing and recombination assignment succeeded, 
   *         false if the map file was malformed, empty, or chromosome mismatch occurred.
   *
   * @note
   * - Assumes `sites` is sorted by position and belongs to a single chromosome.
   * - Sites before the first map entry and after the last entry are extrapolated linearly.
   * - The last site’s recombination probability is always set to zero.
   * - Uses an exponential model for recombination probability conversion:
   *   \f[
   *      \text{recom} = \max\left(\frac{1 - e^{-\Delta / 50}}{2}, \text{recom_min}\right)
   *   \f]
   *   where \f$\Delta\f$ is the genetic distance between consecutive variants.
   *
   * @warning
   * Assertions are used to check for invalid map entries (duplicate positions, overflow).
   * These should be handled gracefully in production code.
   */
  static bool parse_map_file(const std::string& map_file_path, std::vector<target_variant>& sites, float recom_min);

  /**
   * @brief Haldane’s mapping function.
   *
   * Converts a genetic distance (in centiMorgans) into the probability of 
   * recombination between two loci.
   *
   * @param cm Genetic distance in centiMorgans.
   * @return Recombination probability.
   *
   * @f[
   *   r = \frac{1 - e^{-d/50}}{2}
   * @f]
   * where \f$d\f$ is the distance in cM.
   */
  static double haldane(double cm) { return (1. - std::exp(-cm/50.))/2.;}

  /**
   * @brief Convert genetic distance to switch probability (default decay rate).
   *
   * @param cm Genetic distance in centiMorgans.
   * @return Switch probability.
   *
   * @f[
   *   p = 1 - e^{-d/100}
   * @f]
   */
  static double cm_to_switch_prob(double cm) { return 1. - std::exp(-cm/100.);}

  /**
   * @brief Convert genetic distance to switch probability with a custom decay rate.
   *
   * @param cm Genetic distance in centiMorgans.
   * @param decay_rate Custom decay factor scaling the conversion.
   * @return Switch probability.
   *
   * @f[
   *   p = 1 - e^{-\lambda d/100}
   * @f]
   * where \f$\lambda\f$ is the decay rate.
   */
  static double cm_to_switch_prob(double cm, double decay_rate) { return 1. - std::exp(-decay_rate * cm/100.);}

  /**
   * @brief Inverse Haldane’s mapping function.
   *
   * Converts a recombination probability back into genetic distance (in cM).
   *
   * @param recom_prob Recombination probability.
   * @return Genetic distance (cM).
   *
   * @f[
   *   d = 50 \cdot \ln\left(\frac{1}{1 - 2r}\right)
   * @f]
   * where \f$r\f$ is the recombination probability.
   */
  static double haldane_inverse(double recom_prob) { return 50. * std::log(1. / (1. - 2. * recom_prob)); }

  /**
 * @brief Convert switch probability to genetic distance.
 *
 * @param recom_prob Switch probability.
 * @return Genetic distance (cM).
 *
 * @f[
 *   d = 100 \cdot \ln\left(\frac{1}{1 - p}\right)
 * @f]
 * where \f$p\f$ is the switch probability.
 */
  static double switch_prob_to_cm(double recom_prob) { return 100. * std::log(1. / (1. - recom_prob)); }

private:
  /**
   * @brief Read a single entry from a genetic map file.
   *
   * Parses one line from the genetic map input stream into a structured
   * record (`map_file_line`). The function supports both the new format
   * (3 tab-separated columns: chrom, pos, map_value) and the old format
   * (4 columns: chrom, [discard], map_value, pos).
   *
   * @param ifs Input stream of the map file.
   * @param entry Output parameter to store the parsed line (chromosome, position, genetic map value).
   * @param new_format Whether the map file follows the new format (true) or old format (false).
   *
   * @return True if an entry was successfully read, false otherwise.
   *
   * @note An empty chromosome field or failed stream read will cause the function to return false.
   */
  static bool read_entry(std::istream& ifs, map_file_line& entry, bool new_format);
};

/**
 * @brief A reader and interpolator for genetic map files.
 *
 * The `genetic_map_file` class provides an interface to read recombination
 * rate data from a genetic map file and interpolate genetic distances
 * (in centimorgans) for arbitrary variant positions.
 *
 * The file may be in one of two formats:
 * - **New format**: Three columns — `chrom pos cM`
 * - **Legacy format**: Four columns — `chrom <discard> cM pos`
 *
 * Only entries corresponding to the specified target chromosome are processed.
 *
 * @see record, read_record(), interpolate_centimorgan()
 */
class genetic_map_file
{
public:
  /**
   * @brief A single line entry from the genetic map file.
   *
   * Represents a mapping of a genomic coordinate to its genetic distance.
   */
  struct record
  {
    std::string chrom;   ///< Chromosome identifier (e.g., "1", "chrX").
    std::size_t pos = 0; ///< Basepair position on the chromosome.
    double map_value = 0.; ///< Genetic distance (in centimorgans).
  };
private:
  shrinkwrap::istream ifs_;  ///< Input stream for the genetic map file.
  std::string target_chrom_; ///< Chromosome of interest.
  record prev_rec_;          ///< Previously read record for interpolation.
  record cur_rec_;           ///< Current record used for interpolation.
  bool good_;                ///< Status flag: true if file is valid and usable.
  bool new_format_;          ///< True if file follows the new three-column format.
public:
  /**
   * @brief Constructs a genetic_map_file object for reading recombination map data.
   *
   * This constructor opens a genetic map file and prepares it for iterating over
   * recombination map records associated with a specific chromosome.
   *
   * The constructor:
   * - Opens the map file stream.
   * - Detects whether the file is in "new format" (three-column format with tab separators)
   *   or the older PLINK-style format.
   * - Skips header lines (those beginning with `#`).
   * - Searches for the first record matching the requested chromosome.
   * - Reads and buffers the first two records (`prev_rec_` and `cur_rec_`) for iteration.
   *
   * @param map_file_path Path to the genetic map file.
   * @param chrom Chromosome identifier to extract records for.
   *
   * @note
   * - If the file contains invalid header formatting, an error is printed and
   *   the object is marked invalid (`good_ = false`).
   * - If no records for the target chromosome are found, or only a single record
   *   exists, the object is marked invalid.
   * - The validity of the object should be checked with `good()` before use.
   *
   * Example:
   * @code
   * genetic_map_file gmap("genetic_map_chr1.txt", "1");
   * if (!gmap.good()) {
   *   throw std::runtime_error("Failed to load map file");
   * }
   * @endcode
   */
  genetic_map_file(const std::string& map_file_path, const std::string& chrom);


  /**
   * @brief Check whether the genetic map file was loaded successfully.
   *
   * This function returns the internal state flag that indicates whether the
   * constructor succeeded in opening and preparing the file.
   *
   * @return true if the object is valid and ready for use, false otherwise.
   */
  bool good() const { return good_; }

  /**
   * @brief Implicit conversion to bool for validity checking.
   *
   * Provides a shorthand way to test whether the object is in a usable state.
   * Equivalent to calling good().
   *
   * Example:
   * @code
   * genetic_map_file gmap("map.txt", "1");
   * if (gmap) {
   *   // safe to use
   * } else {
   *   // file not valid
   * }
   * @endcode
   *
   * @return true if the object is valid and ready for use, false otherwise.
   */
  operator bool() const { return good_; }

  /**
   * @brief Interpolate the genetic map position (in centimorgans) for a variant.
   *
   * Given a genomic coordinate (basepair position), this method estimates
   * the corresponding centimorgan (cM) value using linear interpolation
   * between records in the genetic map file. The interpolation relies on
   * two adjacent map records (`prev_rec_` and `cur_rec_`).
   *
   * - If the position is before the first record, interpolation assumes
   *   a proportional relationship between basepairs and cM using the
   *   first record as reference.
   * - If the position falls between two known records, interpolation is
   *   performed linearly between their map values.
   * - If the position extends beyond the current block, additional records
   *   are read until the correct interval is found or until the chromosome
   *   ends.
   *
   * @param variant_pos The genomic coordinate (basepair position) of the variant.
   * @return The interpolated centimorgan (cM) value. If the object is invalid
   *         (`good_ == false`), returns NaN.
   *
   * @note Assertions ensure that consecutive map records have different
   *       positions. These cases should ideally be handled more gracefully.
   *
   * @warning Extrapolation beyond the last record assumes constant recombination
   *          rate (`basepair_cm`). Interpret with caution.
   *
   * @see good()
   */
  double interpolate_centimorgan(std::size_t variant_pos);
private:

  /**
   * @brief Read the next record from the genetic map file.
   *
   * This method extracts a single record from the input stream `ifs_` into
   * the provided `entry`. The parsing format depends on whether the map file
   * is in the "new format" (three-column: `chrom pos cM`) or the legacy format
   * (four-column: `chrom <discard> cM pos`).
   *
   * @param[out] entry A reference to a record object that will be populated
   *                   with chromosome, position, and map value.
   *
   * @return `true` if a record was successfully read and is valid,
   *         `false` if the stream is exhausted, corrupted, or the
   *         chromosome field is empty.
   *
   * @note The `discard` column in legacy format is read but ignored.
   * @see interpolate_centimorgan(), good()
   */
  bool read_record(record& rec);
};

#endif // MINIMAC4_RECOMBINATION_HPP