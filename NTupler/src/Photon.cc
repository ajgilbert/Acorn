#include "../interface/Photon.h"

namespace ac {
Photon::Photon()
    : ac::Candidate(),
      isLooseIdPhoton_(false),
      isMediumIdPhoton_(false),
      isTightIdPhoton_(false),
      passElectronVeto_(false),
      hasPixelSeed_(false),
      scEta_(0.),
      hadTowOverEm_(0.),
      full5x5SigmaIetaIeta_(0.),
      chargedIso_(0.),
      neutralHadronIso_(0.),
      photonIso_(0.) {}

Photon::~Photon() {}

void Photon::Print() const { Candidate::Print(); }
}
