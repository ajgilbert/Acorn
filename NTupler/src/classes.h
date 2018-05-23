#include <vector>
#include <map>
#include <utility>
#include <string>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/PileupInfo.h"

namespace { struct dictionary {
  ac::Candidate dummy1;
  std::vector<ac::Candidate> dummy2;
  ac::Muon dummy3;
  std::vector<ac::Muon> dummy4;
  ac::Photon dummy5;
  std::vector<ac::Photon> dummy6;
  ac::GenParticle dummy7;
  std::vector<ac::GenParticle> dummy8;
  std::map<std::string, std::pair<bool, double>> dummy9;
  ac::EventInfo dummy10;
  ac::PileupInfo dummy11;
  std::vector<ac::PileupInfo> dummy12;
};
}

