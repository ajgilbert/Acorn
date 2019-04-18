#include "../interface/PFJet.h"

namespace ac {

PFJet::PFJet()
    : Candidate::Candidate(),
      passesJetID_(false) {}

  PFJet::~PFJet() {}

  void PFJet::Print() const {
    Candidate::Print();
  }
}
