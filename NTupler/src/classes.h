#include <vector>
#include <map>
#include <utility>
#include <string>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/EventInfo.h"

namespace { struct dictionary {
  ac::Candidate dummy1;
  std::vector<ac::Candidate> dummy2;
  ac::GenParticle dummy8;
  std::vector<ac::GenParticle> dummy9;
  std::map<std::string, std::pair<bool, double>> dummy3;
  ac::EventInfo dummy10;

};
}

