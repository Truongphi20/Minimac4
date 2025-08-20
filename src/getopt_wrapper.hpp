#ifndef GETOPT_WRAPPER_HPP
#define GETOPT_WRAPPER_HPP

#include <cstdint>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <getopt.h>

/**
 * @class getopt_wrapper
 * @brief A wrapper for parsing command-line options using getopt and getopt_long.
 *
 * This class simplifies handling both short and long command-line options,
 * automatically generating the short option string and formatting help output.
 *
 * @details
 * The class provides:
 * - `option_with_desc`: A structure extending `struct option` to include
 *   a description string for help messages.
 * - Storage for long options (`long_options_`) and original option
 *   definitions (`opts_`).
 * - The program usage string (`usage_str_`) and the generated short option
 *   string (`short_opt_string_`).
 * - Tracking of the maximum long option length (`max_long_opt_length_`)
 *   for nicely aligned help output.
 */
class getopt_wrapper
{
public:
  /**
   * @struct getopt_wrapper::option_with_desc
   * @brief Extension of `struct option` that includes a description string.
   *
   * This structure inherits from the standard GNU `struct option`
   * (used with `getopt_long`) and adds an additional field for
   * a human-readable description. It is useful for generating
   * formatted help/usage messages alongside command-line parsing.
   *
   * @details
   * Each option includes:
   * - `name`: The long option name (e.g., `"--input"`).
   * - `has_arg`: Whether the option requires an argument 
   *   (`no_argument`, `required_argument`, or `optional_argument`).
   * - `flag`: Pointer to an int variable that is set if the option is found, 
   *   or `nullptr` if the option’s value should be returned.
   * - `val`: The short option character equivalent (e.g., `'i'`).
   * - `description`: A string describing the option, shown in help output.
   *
   * Example:
   * @code
   * getopt_wrapper::option_with_desc opt("input", required_argument, nullptr, 'i',
   *                                      "Path to input file");
   * @endcode
   */
  struct option_with_desc : public ::option
  {
    /**
     * @brief Description of the command-line option for help output.
     */
    const char* description;

    /**
     * @brief Construct an option with description.
     *
     * @param _name        Long option name (without leading dashes).
     * @param _has_arg     Argument requirement (no_argument, required_argument, optional_argument).
     * @param _flag        Pointer to variable set if option is found, or nullptr.
     * @param _val         Value returned (or short option character).
     * @param _description Human-readable description of the option.
     */
    option_with_desc(const char* _name, int _has_arg, int* _flag, int _val, const char* _description)
    {
      name = _name;
      has_arg = _has_arg;
      flag = _flag;
      val = _val;
      description = _description;
    }
  };
protected:
  /**
   * @brief Storage of long options for getopt parsing.
   *
   * This vector holds the `struct option` objects used by `getopt_long`.
   * It is automatically populated based on the provided
   * `option_with_desc` definitions.
   */
  std::vector<option> long_options_;

  /**
   * @brief Original option definitions with descriptions.
   *
   * This vector stores the extended option structures
   * (`option_with_desc`), which include descriptions for generating
   * nicely formatted help/usage messages.
   */
  std::vector<option_with_desc> opts_;

  /**
   * @brief Program usage/help string.
   *
   * Contains a short description of the program usage,
   * typically displayed at the top of the `--help` output.
   */
  std::string usage_str_;

  /**
   * @brief Concatenated short option string.
   *
   * Generated automatically from the provided options.
   * This string is passed to `getopt()` for parsing
   * short options (e.g., `"hv:o:"`).
   */
  std::string short_opt_string_;

  /**
   * @brief Maximum long option length.
   *
   * Used when formatting help output so that all descriptions
   * are aligned to the same column, regardless of option name length.
   */
  std::size_t max_long_opt_length_ = 0;
public:
  /**
   * @brief Constructs a getopt_wrapper object for parsing command-line options.
   *
   * This constructor initializes both short and long options for use with
   * `getopt` and `getopt_long`. It sets up internal data structures for
   * option parsing and generates the short option string automatically.
   *
   * @param usage_str A string describing the program usage, typically
   *        displayed when the user requests help.
   * @param long_opts A vector of `option_with_desc` structures describing
   *        long options, their short equivalents, argument requirements,
   *        and descriptions.
   *
   * @details
   * The constructor performs the following steps:
   * 1. Moves the usage string and option vector into member variables.
   * 2. Resizes the internal `long_options_` array to accommodate all options
   *    plus a terminating zeroed element (required by `getopt_long`).
   * 3. Copies each option from `opts_` into `long_options_`.
   * 4. Updates `max_long_opt_length_` to track the length of the longest
   *    long option name (used for formatting help output).
   * 5. Builds the short option string `short_opt_string_`:
   *    - Adds the short option character if `val` is set.
   *    - Appends `:` if the option requires an argument (`required_argument`).
   */
  getopt_wrapper(std::string usage_str, std::vector<option_with_desc> long_opts) :
    usage_str_(std::move(usage_str)),
    opts_(std::move(long_opts))
  {

    long_options_.resize(opts_.size() + 1, {0, 0, 0, 0});
    auto lit = long_options_.begin();
    for (auto it = opts_.begin(); it != opts_.end(); ++it)
    {
      *(lit++) = *it;
      max_long_opt_length_ = std::max(max_long_opt_length_, it->name ? std::strlen(it->name) : 0);
      if (it->val)
      {
        short_opt_string_ += (char)it->val;
        if (it->has_arg == required_argument)
          short_opt_string_ += ':';
      }
    }
  }

  /**
   * @brief Prints the program usage and all available options to a stream.
   *
   * This function formats and displays the usage string followed by a
   * nicely aligned list of command-line options and their descriptions.
   * Both short and long options are displayed if available.
   *
   * @param os The output stream to which the usage information is written
   *           (e.g., `std::cout` or `std::cerr`).
   *
   * @details
   * The method performs the following steps:
   * 1. Prints the main usage string (`usage_str_`) followed by a blank line.
   * 2. Iterates over all options in `opts_`.
   * 3. For each option with a description:
   *    - Prints the short option (if printable) prefixed with `-`.
   *    - Prints the long option name (if present) prefixed with `--`.
   *    - Aligns descriptions based on the longest long option (`max_long_opt_length_`)
   *      for a clean, columnar format.
   * 4. Flushes the output stream to ensure all text is written.
   *
   * Example output:
   * @code
   * Usage: myprogram [options]
   * 
   *  -h, --help        Display this help message
   *  -v, --version     Show version information
   * @endcode
   */
  void print_usage(std::ostream& os)
  {
    os << usage_str_ << '\n';
    os << '\n';
    for (auto it = opts_.begin(); it != opts_.end(); ++it)
    {
      if (!it->description)
        continue;

      if (std::isprint(it->val))
      {
        if (it->name)
          os << " -" << (char)it->val << ", ";
        else
          os << " -" << (char)it->val << "  ";
      }
      else
        os << "     ";

      std::size_t n_spaces = 2;
      if (it->name)
        os << "--" << it->name;
      else
        n_spaces += 2;

      n_spaces += max_long_opt_length_ - std::strlen(it->name);
      for (std::size_t i = 0; i < n_spaces; ++i)
        os.put(' ');
      os << it->description << '\n';
    }

    os << std::flush;
  }
};

#endif // GETOPT_WRAPPER_HPP
