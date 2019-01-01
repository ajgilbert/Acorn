#include "../interface/Met.h"
#include <map>
#include <string>
#include <vector>

namespace ac {
Met::Met() : Candidate::Candidate(), level_(-1), shift_(-1), sumEt_(0.) {}

Met::~Met() {}

void Met::Print() const { Candidate::Print(); }

}  // namespace ac
