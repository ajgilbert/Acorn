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
  float phi1_;
  float lhe_phi1_;
  unsigned type_;
  float weight_C3w_0p1_;
  float weight_C3w_0p2_;
  float weight_C3w_0p4_;
  float weight_C3w_0p8_;
  float weight_C3w_1p6_;
  float w_pt_;
  float g_pt_;
  float g_eta_;
  float n_pt_;
  float n_eta_;
  float l_pt_;
  float l_eta_;
  int l_charge_;
  float l_g_dr_;
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

  bool IsChargedLepton(ac::GenParticle const& p) const;
  bool IsNeutrino(ac::GenParticle const& p) const;
  bool IsLepton(ac::GenParticle const& p) const;
  bool IsPhoton(ac::GenParticle const& p) const;
};
}

#endif
