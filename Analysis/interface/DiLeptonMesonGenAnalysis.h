#ifndef Acorn_Analysis_DiLeptonMesonGenAnalysis_h
#define Acorn_Analysis_DiLeptonMesonGenAnalysis_h
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

class DiLeptonMesonGenAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(DiLeptonMesonGenAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(DiLeptonMesonGenAnalysis, unsigned, year)

  TTree* tree_;
  float pt_1_;
  float pt_2_;
  float eta_1_;
  float eta_2_;
  float m_ll_;
  float pt_ll_;
  float dr_ll_;
  float gen_trk_pt_1_;
  float gen_trk_pt_2_;
  float gen_trk_eta_1_;
  float gen_trk_eta_2_;
  float gen_trk_phi_1_;
  float gen_trk_phi_2_;
  float reco_trk_pt_1_;
  float reco_trk_pt_2_;
  float reco_trk_eta_1_;
  float reco_trk_eta_2_;
  float reco_trk_phi_1_;
  float reco_trk_phi_2_;
  float reco_trk_charge_1_;
  float reco_trk_charge_2_;
  float reco_trk_dR_;
  float reco_trk_mass_;
  float reco_trk_pt_;
  float reco_trk_iso_;
  float reco_trk_looser_iso_;
  unsigned reco_trk_id_1_;
  unsigned reco_trk_id_2_;
  float reco_higgs_mass_;
  float reco_higgs_pt_;

 public:
  DiLeptonMesonGenAnalysis(std::string const& name);
  virtual ~DiLeptonMesonGenAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
