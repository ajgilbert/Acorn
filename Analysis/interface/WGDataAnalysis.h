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
  CLASS_MEMBER(WGDataAnalysis, std::string, gen_classify)

  LookupFilter filters_IsoMu24_;
  LookupFilter filters_IsoTkMu24_;
  LookupFilter filters_IsoMu27_;

  std::shared_ptr<RooWorkspace> ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;

  TTree* tree_;

  // truth properties
  unsigned gen_proc_;

  unsigned n_m_; // number of medium muons

  // m0: Main muon variables
  float m0_pt_;
  float m0_eta_;
  float m0_phi_;
  float m0_iso_;
  bool m0_trg_; // trigger fired and object matched

  // m1: Second muon variables
  float m1_pt_;
  float m1_eta_;
  float m1_phi_;
  float m1_iso_;

  // di-muon variables
  float m0m1_M_;
  float m0m1_dr_;
  bool m0m1_os_;

  // number of reco'd photons (no ID/Iso beyond miniaod presel)
  unsigned n_p_;

  // p0: Main photon variables, defined if n_p >= 1
  float p0_pt_;
  float p0_eta_;
  float p0_phi_;
  float p0_chiso_;
  float p0_neiso_;
  float p0_phiso_;
  float p0_hovere_;
  float p0_sigma_;
  bool p0_haspix_;
  bool p0_medium_noch_;
  bool p0_medium_; // also passes medium ID
  bool p0_tight_; // also passes tight ID
  bool p0_isprompt_; // truth matched to prompt photon

  // met
  float met_;
  float met_phi_;

  // composite variables
  float m0met_mt_;

  // composite variables defined inf n_p >= 1
  float m0p0_dr_;
  float m0p0_dphi_;
  float m0p0_M_;
  float reco_phi_;

  // vetos
  unsigned n_vm_; // number of additional veto muons
  float vm_p0_dr_; // if n_vm >= 1, deltaR between photon and closest veto muon

  // event weights
  float wt_def_; // default weight
  float wt_pu_;
  float wt_m0_; // trk/ID/Iso weight for m0
  float wt_trg_m0_; // trigger weight for m0
  float wt_m1_; // trk/ID/Iso weight for m1
  float wt_p0_; // ID/iso weight for p0

  mutable TRandom3 rng;


 public:
  WGDataAnalysis(std::string const& name);
  virtual ~WGDataAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

  void SetDefaults();
};
}  // namespace ac

#endif
