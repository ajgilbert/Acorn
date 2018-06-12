#include "../interface/Met.h"
#include <map>
#include <string>
#include <vector>

namespace ac {
Met::Met() : Candidate::Candidate(), sumEt_(0.) {}

Met::~Met() {}

void Met::Print() const { Candidate::Print(); }

}  // namespace ac
