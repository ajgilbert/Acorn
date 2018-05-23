#include "../interface/PileupInfo.h"

namespace ac {
PileupInfo::PileupInfo() : numInteractions_(0), bunchCrossing_(0), trueNumInteractions_(0.) {}

PileupInfo::~PileupInfo() {}

void PileupInfo::Print() const {}
}  // namespace ac
