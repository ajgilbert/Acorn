#ifndef Acorn_Analysis_DiMuonAnalysis_h
#define Acorn_Analysis_DiMuonAnalysis_h
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

class DiMuonAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(DiMuonAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(DiMuonAnalysis, unsigned, year)
  CLASS_MEMBER(DiMuonAnalysis, bool, is_data)
  CLASS_MEMBER(DiMuonAnalysis, std::string, corrections)

  LookupFilter filters_IsoMu24_;
  LookupFilter filters_IsoTkMu24_;
  LookupFilter filters_IsoMu27_;

  std::shared_ptr<RooWorkspace> ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;

  TTree* tree_;
  float pt_1_;
  float pt_2_;
  float eta_1_;
  float eta_2_;
  float m_ll_;
  float pt_ll_;
  float dr_ll_;
  bool trg_1_;
  bool trg_2_;
  float wt_pu_;
  float wt_1_;
  float wt_2_;
  float wt_trg1_;
  float wt_trg2_;

 public:
  DiMuonAnalysis(std::string const& name);
  virtual ~DiMuonAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
