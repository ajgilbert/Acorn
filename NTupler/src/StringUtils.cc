#include "Acorn/NTupler/interface/StringUtils.h"

#include <string>
#include <vector>
#include "boost/algorithm/string.hpp"

namespace ac {

std::string TrimString(std::string const& input) {
  std::string input_copy = input;
  boost::trim(input_copy);
  return input_copy;
}


std::vector<std::string> TrimAndSplitString(std::string const& input) {
  std::string input_copy = input;
  boost::trim(input_copy);
  std::vector<std::string> split;
  boost::split(split, input_copy, boost::is_any_of(" \t"), boost::token_compress_on);
  return split;
}
}  // namespace ac
