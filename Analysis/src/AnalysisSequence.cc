#include "Acorn/Analysis/interface/AnalysisSequence.h"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include "Acorn/Analysis/interface/ModuleBase.h"

namespace ac {

void Sequence::InsertSequence(std::string name, ac::AnalysisBase& ana) {
  for (auto m : seq) {
    ana.AddModule(name, m.get());
  }
}

}
