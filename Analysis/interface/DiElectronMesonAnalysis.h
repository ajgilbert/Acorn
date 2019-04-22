#ifndef Acorn_Analysis_DiElectronMesonAnalysis_h
#define Acorn_Analysis_DiElectronMesonAnalysis_h
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

class DiElectronMesonAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(DiElectronMesonAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(DiElectronMesonAnalysis, unsigned, year)
  CLASS_MEMBER(DiElectronMesonAnalysis, bool, is_data)
  CLASS_MEMBER(DiElectronMesonAnalysis, std::string, corrections)

  LookupFilter filters_Ele35_;
  LookupFilter filters_Ele32_;
  LookupFilter filters_Ele27_;

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
  float wt_trg_;
  float eff_trg1_data_;
  float eff_trg2_data_;
  float eff_trg1_mc_;
  float eff_trg2_mc_;
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
  DiElectronMesonAnalysis(std::string const& name);
  virtual ~DiElectronMesonAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
