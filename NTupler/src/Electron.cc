#include "../interface/Electron.h"

namespace ac {

Electron::Electron()
    : Candidate::Candidate(),
      isCutBasedVetoElectron_(false),
      isCutBasedLooseElectron_(false),
      isCutBasedMediumElectron_(false),
      isCutBasedTightElectron_(false),
      isMVAwp90Electron_(false),
      isMVAwp80Electron_(false),
      isHEEPElectron_(false),
      dxy_(0.),
      dz_(0.),
      scEta_(0.),
      scEnergy_(0.){}

  Electron::~Electron() {}

  void Electron::Print() const {
    Candidate::Print();
  }
}
