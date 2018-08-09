#ifndef Acorn_Analysis_WGDataAnalysis_h
#define Acorn_Analysis_WGDataAnalysis_h
#include <string>
#include <cstdint>
#include "boost/range/algorithm_ext/erase.hpp"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include "TRandom3.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/GenParticle.h"

namespace ac {

class WGDataAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(WGDataAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(WGDataAnalysis, unsigned, year)
  CLASS_MEMBER(WGDataAnalysis, bool, is_data)
  CLASS_MEMBER(WGDataAnalysis, std::string, corrections)

  LookupFilter filters_IsoMu24_;
  LookupFilter filters_IsoTkMu24_;
  LookupFilter filters_IsoMu27_;

  std::shared_ptr<RooWorkspace> ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;

  TTree* tree_;
  float pt_m_;
  float pt_p_;
  float eta_m_;
  float eta_p_;
  float pt_met_;
  float mt_;
  bool trg_m_;
  float wt_pu_;
  float wt_m_;
  float wt_trg_m_;

 public:
  WGDataAnalysis(std::string const& name);
  virtual ~WGDataAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
