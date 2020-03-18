#include "../interface/Muon.h"
// #include "../interface/city.h"

namespace ac {

Muon::Muon()
    : Candidate::Candidate(),
      isLooseMuon_(false),
      isMediumMuon_(false),
      isTightMuon_(false),
      dxy_(0.),
      dz_(0.),
      pfIsoSumChargedHadronPt_(0.),
      pfIsoSumNeutralHadronEt_(0.),
      pfIsoSumPhotonEt_(0.),
      pfIsoSumPUPt_(0.) {}

  Muon::~Muon() {}

  void Muon::Print() const {
    Candidate::Print();
  }
}
