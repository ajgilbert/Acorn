#include "../interface/Photon.h"

namespace ac {
Photon::Photon()
    : ac::Candidate(),
      isLooseIdPhoton_(false),
      isMediumIdPhoton_(false),
      isTightIdPhoton_(false),
      passElectronVeto_(false),
      hasPixelSeed_(false) {}

Photon::~Photon() {}

void Photon::Print() const { Candidate::Print(); }
}
