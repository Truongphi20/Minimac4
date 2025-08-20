#ifndef MINIMAC4_VARIANT_HPP
#define MINIMAC4_VARIANT_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <limits>

/**
 * @struct target_variant
 * @brief Represents a variant in the target dataset.
 *
 * Contains information about the variant's chromosome, position, alleles,
 * genotype data, and flags indicating whether the site exists in the target
 * or reference file.
 */
struct target_variant
{
  std::string chrom;      ///< Chromosome name
  std::uint32_t pos;      ///< 1-based position on the chromosome
  std::string id;         ///< Variant identifier
  std::string ref;        ///< Reference allele
  std::string alt;        ///< Alternate allele
  bool in_tar;            ///< True if site exists in the target file
  bool in_ref;            ///< True if site exists in the reference file
  float af;               ///< Allele frequency
  float err;              ///< Error rate
  float recom;            ///< Recombination rate
  std::vector<std::int8_t> gt; ///< Genotype data for each sample
};

/**
 * @struct reference_site_info
 * @brief Stores information about a site in the reference dataset.
 *
 * Includes basic variant information as well as error and recombination rates,
 * and centimorgan position if available.
 */
struct reference_site_info
{
  std::string chrom;      ///< Chromosome name
  std::uint32_t pos = 0;  ///< Position on the chromosome
  std::string id;         ///< Variant identifier
  std::string ref;        ///< Reference allele
  std::string alt;        ///< Alternate allele
  float err = std::numeric_limits<float>::quiet_NaN(); ///< Error rate
  float recom = std::numeric_limits<float>::quiet_NaN(); ///< Recombination rate
  double cm = std::numeric_limits<double>::quiet_NaN(); ///< Centimorgan position

  reference_site_info() {}

  reference_site_info(std::string _chrom,
                      std::uint32_t _pos,
                      std::string _id,
                      std::string _ref,
                      std::string _alt,
                      float _err,
                      float _recom,
                      double _cm)
    : chrom(std::move(_chrom)),
      pos(_pos),
      id(_id),
      ref(std::move(_ref)),
      alt(std::move(_alt)),
      err(_err),
      recom(_recom),
      cm(_cm)
  {}
};

/**
 * @struct reference_variant
 * @brief Extends reference_site_info to include genotype data and allele counts.
 */
struct reference_variant : public reference_site_info
{
  reference_variant() {}

  reference_variant(const std::string& _chrom,
                    std::uint32_t _pos,
                    const std::string& _id,
                    const std::string& _ref,
                    const std::string& _alt,
                    float _err,
                    float _recom,
                    double _cm,
                    std::size_t _ac,
                    std::vector<std::int8_t> _gt)
    : reference_site_info(_chrom, _pos, _id, _ref, _alt, _err, _recom, _cm),
      ac(_ac),
      gt(std::move(_gt))
  {}

  std::size_t ac = 0;                   ///< Allele count
  std::vector<std::int8_t> gt;          ///< Genotype vector
};

/**
 * @struct sparse_ref_variant
 * @brief Represents a reference variant with sparse allele information.
 *
 * Stores allele count and offsets to alternate alleles in addition to the basic
 * reference site info.
 */
struct sparse_ref_variant : public reference_site_info
{
  std::size_t ac;                        ///< Allele count
  std::vector<std::uint32_t> alt_allele_offsets; ///< Offsets of alternate alleles

  sparse_ref_variant(const std::string& _chrom,
                     std::uint32_t _pos,
                     const std::string& _id,
                     const std::string& _ref,
                     const std::string& _alt,
                     float _err,
                     float _recom,
                     double _cm,
                     std::size_t _ac,
                     const std::size_t* off_it, const std::size_t* off_it_end)
    : reference_site_info(_chrom, _pos, _id, _ref, _alt, _err, _recom, _cm),
      ac(_ac),
      alt_allele_offsets(off_it, off_it_end)
  {}
};

#endif // MINIMAC4_VARIANT_HPP
