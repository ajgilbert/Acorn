#include "../interface/PFJet.h"

namespace ac {

PFJet::PFJet()
    : Candidate::Candidate(),
      passesJetID_(false),
      deepCSVDiscriminatorBvsAll_(0.) {}

  PFJet::~PFJet() {}

  void PFJet::Print() const {
    Candidate::Print();
  }
}
