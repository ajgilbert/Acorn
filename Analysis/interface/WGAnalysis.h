#ifndef Acorn_Analysis_WGAnalysis_h
#define Acorn_Analysis_WGAnalysis_h
#include <string>
#include <cstdint>
#include <memory>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "TRandom3.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/Analysis/interface/StandaloneReweight.h"
#include "Acorn/NTupler/interface/GenParticle.h"

namespace ac {

class WGAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(WGAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(WGAnalysis, bool, add_standalone)
  CLASS_MEMBER(WGAnalysis, bool, add_rivet)

  TTree * tree_;

  // Variables in common with RIVET
  bool is_wg_gen_;
  float gen_l0_pt_;
  float gen_l0_eta_;
  int gen_l0_q_;
  unsigned gen_pdgid_;
  float gen_p0_pt_;
  float gen_p0_eta_;
  bool gen_p0_frixione_;
  float gen_l0p0_dr_;
  float gen_n0_pt_;
  float gen_n0_eta_;
  float gen_met_pt_;
  float gen_true_phi_;
  float gen_true_phi_f_;

  bool is_wg_riv_;
  float riv_l0_pt_;
  float riv_l0_eta_;
  int riv_l0_q_;
  unsigned riv_pdgid_;
  float riv_p0_pt_;
  float riv_p0_eta_;
  bool riv_p0_frixione_;
  float riv_l0p0_dr_;
  float riv_n0_pt_;
  float riv_n0_eta_;
  float riv_met_pt_;
  float riv_true_phi_;
  float riv_true_phi_f_;


  float gen_phi_;
  float gen_phi_f_;
  float lhe_true_phi_;
  float lhe_true_phi_f_;
  unsigned type_;
  float wt_def_;
  float wt_C3w_0p0_;
  float wt_C3w_0p1_;
  float wt_C3w_0p2_;
  float wt_C3w_0p4_;
  float wt_C3w_0p67_;
  float wt_C3w_1p0_;
  float gen_W_pt_;
  float gen_l0p0_dphi_;
  float gen_Wp0_dphi_;
  float j0_pt_;
  float j0_eta_;
  int n_jets_;
  unsigned nparts_;
  bool valid_mt_;

  // Weights computed using standalone ME RW
  float st_C3w_0p0_;
  float st_C3w_0p1_;
  float st_C3w_0p2_;
  float st_C3w_0p4_;
  float st_C3w_0p67_;
  float st_C3w_1p0_;

  mutable TRandom3 rng;

  std::shared_ptr<StandaloneReweight> rw_;
  std::map<int, int> replace_pdg_;

 public:
  WGAnalysis(std::string const& name);
  virtual ~WGAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
  void SetDefaults();
};
}

#endif
