#ifndef Acorn_Analysis_AnalysisSequence_h
#define Acorn_Analysis_AnalysisSequence_h

#include <vector>
#include <stdexcept>
#include <map>
#include <string>
#include <iostream>
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/Analysis/interface/AnalysisBase.h"

namespace ac {

class Sequence {
 private:
  std::vector<std::shared_ptr<ac::ModuleBase>> seq;

 public:
  typedef std::vector<std::shared_ptr<ac::ModuleBase>> ModuleSequence;
  Sequence() = default;
  ~Sequence() = default;

  template<class T>
  void BuildModule(T const& mod) {
     seq.push_back(std::make_shared<T>(mod));
  }

  void InsertSequence(std::string name, ac::AnalysisBase & ana);
};

}
#endif
