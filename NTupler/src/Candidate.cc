#include "../interface/Candidate.h"

namespace ac {
// Constructors/Destructors
Candidate::Candidate() : vector_(Vector()), id_(0), charge_(0) {}

Candidate::~Candidate() {}

void Candidate::Print() const {
  std::cout << "[pt,eta,phi,M] = " << vector_ << " charge = " << charge_
            << std::endl;
}
}
