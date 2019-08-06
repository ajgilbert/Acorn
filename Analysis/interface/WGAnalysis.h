#ifndef Acorn_Analysis_WGAnalysis_h
#define Acorn_Analysis_WGAnalysis_h
#include <string>
#include <cstdint>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include "TRandom3.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/NTupler/interface/GenParticle.h"

namespace ac {

class WGAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(WGAnalysis, fwlite::TFileService*, fs)

  TTree * tree_;
  float gen_phi_;
  float gen_phi_f_;
  float true_phi_;
  float true_phi_f_;
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
  float w_pt_;
  float g_pt_;
  float g_eta_;
  float n_pt_;
  float n_eta_;
  float l_pt_;
  float l_eta_;
  int l_charge_;
  float l_g_dr_;
  float l_g_dphi_;
  float w_g_dphi_;
  unsigned nparts_;
  bool valid_mt_;
  mutable TRandom3 rng;

 public:
  WGAnalysis(std::string const& name);
  virtual ~WGAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}

#endif
