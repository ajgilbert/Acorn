#ifndef Acorn_NTupler_StringUtils_h
#define Acorn_NTupler_StringUtils_h

#include<vector>
#include<string>

namespace ac {

std::string TrimString(std::string const& input);

std::vector<std::string> TrimAndSplitString(std::string const& input);
}

#endif
