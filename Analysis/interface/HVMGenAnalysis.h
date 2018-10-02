#ifndef Acorn_Analysis_HVMGenAnalysis_h
#define Acorn_Analysis_HVMGenAnalysis_h

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

class HVMGenAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(HVMGenAnalysis, fwlite::TFileService*, fs)

  TTree * tree_;
  float muon1_pt_; 
  float muon2_pt_; 
  float muon1_eta_; 
  float muon2_eta_; 
  float muon1_phi_; 
  float muon2_phi_; 
  float pion1_pt_;
  float pion2_pt_;
  float pion3_pt_;
  float pion4_pt_;
  float pion5_pt_;
  float pion1_eta_;
  float pion2_eta_;
  float pion3_eta_;
  float pion4_eta_;
  float pion5_eta_;
  float pion1_phi_;
  float pion2_phi_;
  float pion3_phi_;
  float pion4_phi_;
  float pion5_phi_;
  float pion1_charge_;
  float pion2_charge_;
  float pion3_charge_;
  float pion4_charge_;
  float pion5_charge_;
  bool pion1_fromres_;
  bool pion2_fromres_;
  bool pion3_fromres_;
  bool pion4_fromres_;
  bool pion5_fromres_;
  float res_mass_;
  float res_pt_;
  float res_eta_;
  float res_phi_;
  float genres_mass_;
  float genres_pt_;
  float genres_eta_;
  float genres_phi_;
  float pions_dr_;


 public:
  HVMGenAnalysis(std::string const& name);
  virtual ~HVMGenAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

};
}

#endif
