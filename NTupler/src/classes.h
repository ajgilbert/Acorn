#include <vector>
#include <map>
#include <utility>
#include <string>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/NTupler/interface/RunLumiMap.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/PFJet.h"

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
  ac::TriggerObject dummy13;
  std::vector<ac::TriggerObject> dummy14;
  RunLumiMap dummy15;
  std::map<unsigned, std::set<unsigned>> dummy16;
  ac::Met dummy17;
  std::vector<ac::Met> dummy18;
  ac::Track dummy19;
  std::vector<ac::Track> dummy20;
  ac::Electron dummy21;
  std::vector<ac::Electron> dummy22;
  ac::PFJet dummy23;
  std::vector<ac::PFJet> dummy24;
};
}

