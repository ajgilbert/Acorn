#ifndef Acorn_Analysis_HVMTagAndProbe_h
#define Acorn_Analysis_HVMTagAndProbe_h
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

class HVMTagAndProbe : public ModuleBase {
 private:
  CLASS_MEMBER(HVMTagAndProbe, fwlite::TFileService*, fs)
  CLASS_MEMBER(HVMTagAndProbe, unsigned, year)
  CLASS_MEMBER(HVMTagAndProbe, bool, is_data)
  CLASS_MEMBER(HVMTagAndProbe, std::string, corrections)

  LookupFilter filters_IsoMu24_;
  LookupFilter filters_IsoTkMu24_;
  LookupFilter filters_IsoMu27_;

  std::shared_ptr<RooWorkspace> ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;

  TTree* tree_;

  unsigned run_;

  unsigned n_vtx_; // number of reco vertices

  float t_pt_;
  float t_eta_;
  float t_phi_;
  int t_q_;
  bool t_id_;
  bool t_rand_;

  float p_pt_;
  float p_eta_;
  float p_phi_;
  bool p_id_;
  int p_q_;

  float m_ll_;

  bool t_trg_;
  bool p_trg_;

  float p_trk_iso_;

  // event weights
  float wt_def_; // default weight

  mutable TRandom3 rng;

 public:
  HVMTagAndProbe(std::string const& name);
  virtual ~HVMTagAndProbe();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

  void SetDefaults();

  bool PassesTrigger(ac::Muon const* e, ac::TreeEvent * event) const;
};
}  // namespace ac

#endif
