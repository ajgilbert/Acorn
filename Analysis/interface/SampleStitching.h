#ifndef Acorn_Analysis_SampleStitching_h
#define Acorn_Analysis_SampleStitching_h
#include <string>
#include <cstdint>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include "TRandom3.h"
#include "Acorn/NTupler/interface/json.hpp"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/NTupler/interface/GenParticle.h"

namespace ac {

class SampleStitching : public ModuleBase {
 public:
  struct SampleInfo {
    std::string name;
    bool target;
    std::vector<bool> has_min;
    std::vector<bool> has_max;
    std::vector<double> min;
    std::vector<double> max;
    unsigned events;
    double xsec;
  };

 private:
  CLASS_MEMBER(SampleStitching, fwlite::TFileService*, fs)

  TTree* tree_;
  std::vector<double> vars_;
  std::vector<std::string> binned_;
  std::vector<SampleInfo> samples_;

 public:
  SampleStitching(std::string const& name, nlohmann::json const& cfg);
  virtual ~SampleStitching();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
