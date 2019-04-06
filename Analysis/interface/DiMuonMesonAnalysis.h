#ifndef Acorn_Analysis_DiMuonMesonAnalysis_h
#define Acorn_Analysis_DiMuonMesonAnalysis_h
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

class DiMuonMesonAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(DiMuonMesonAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(DiMuonMesonAnalysis, unsigned, year)
  CLASS_MEMBER(DiMuonMesonAnalysis, bool, is_data)
  CLASS_MEMBER(DiMuonMesonAnalysis, std::string, corrections)

  LookupFilter filters_IsoMu24_;
  LookupFilter filters_IsoTkMu24_;
  LookupFilter filters_IsoMu27_;

  std::shared_ptr<RooWorkspace> ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;

  TTree* tree_;
  int  electron_overlaps_;
  int  muon_overlaps_;
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
  float wt_rhoiso_;
  float highestpt_pair_iso_;
  float highestpt_pair_looser_iso_;
  float highestpt_pair_dR_;
  float highestpt_pair_mass_;
  float highestpt_pair_pt_;
  float highestpt_pair_eta_;
  float highestpt_pair_phi_;
  float highestpt_pair_1_phi_;
  float highestpt_pair_1_eta_;
  float highestpt_pair_1_pt_;
  float highestpt_pair_2_phi_;
  float highestpt_pair_2_eta_;
  float highestpt_pair_2_pt_;
  unsigned highestpt_pair_id_1_;
  unsigned highestpt_pair_id_2_;
  float highestpt_pair_reco_higgs_mass_;
  float highestpt_pair_reco_higgs_pt_;
  int nAddEle_;
  float Zrho_dphi_;

 public:
  DiMuonMesonAnalysis(std::string const& name);
  virtual ~DiMuonMesonAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
