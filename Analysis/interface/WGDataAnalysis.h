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
#include "Acorn/Analysis/interface/RoccoR.h"

namespace ac {

class WGDataAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(WGDataAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(WGDataAnalysis, unsigned, year)
  CLASS_MEMBER(WGDataAnalysis, bool, is_data)
  CLASS_MEMBER(WGDataAnalysis, std::string, corrections)
  CLASS_MEMBER(WGDataAnalysis, std::string, gen_classify)
  CLASS_MEMBER(WGDataAnalysis, bool, do_wg_gen_vars)
  CLASS_MEMBER(WGDataAnalysis, bool, check_is_zg)
  CLASS_MEMBER(WGDataAnalysis, bool, check_is_wwg)
  CLASS_MEMBER(WGDataAnalysis, bool, do_presel)
  CLASS_MEMBER(WGDataAnalysis, int, var_set)
  CLASS_MEMBER(WGDataAnalysis, int, correct_e_energy)
  CLASS_MEMBER(WGDataAnalysis, int, correct_p_energy)
  CLASS_MEMBER(WGDataAnalysis, int, correct_m_energy)
  CLASS_MEMBER(WGDataAnalysis, int, shift_met)
  CLASS_MEMBER(WGDataAnalysis, int, scale_weights)
  CLASS_MEMBER(WGDataAnalysis, int, pdf_begin)
  CLASS_MEMBER(WGDataAnalysis, int, pdf_end)
  CLASS_MEMBER(WGDataAnalysis, std::string, rc_file)

  LookupFilter filters_IsoMu24_;
  LookupFilter filters_IsoTkMu24_;
  LookupFilter filters_IsoMu27_;
  LookupFilter filters_Mu50_;

  LookupFilter filters_Ele27_;
  LookupFilter filters_Ele32_;
  LookupFilter filters_Ele32_L1DoubleEG_;
  LookupFilter filters_Ele32_L1DoubleEG_seed_;

  std::shared_ptr<RooWorkspace> ws_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;

  RoccoR rc_;

  TTree* tree_;

  unsigned run_;

  // truth properties
  unsigned gen_proc_;
  bool gen_is_zg_; // Event is in the phase space covered by the NLO ZG sample

  unsigned n_vtx_; // number of reco vertices
  unsigned metfilters_; // metfilter bits (0 == OK)

  unsigned n_pre_m_; // number of medium/tight muons
  unsigned n_pre_e_; // number of medium/tight electrons
  unsigned n_veto_m_; // number of veto muons
  unsigned n_veto_e_; // number of veto electrons

  // l0: Main muon/electron variables
  float l0_pt_;
  float l0_eta_;
  float l0_phi_;
  float l0_iso_;
  float l0_tkiso_;
  bool l0_tight_;
  bool l0_trg_; // trigger fired and object matched
  bool l0_trg_2_; // 2nd trigger option
  int l0_q_;
  unsigned l0_pdgid_;

  // m1: Second muon variables
  float l1_pt_;
  float l1_eta_;
  float l1_phi_;
  float l1_iso_;

  // di-muon variables
  float l0l1_M_;
  float l0l1_dr_;
  float l0l1_pt_;
  bool l0l1_os_;

  // number of reco'd photons (no ID/Iso beyond miniaod presel)
  unsigned n_pre_p_;

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
  bool p0_eveto_;
  bool p0_loose_;
  bool p0_medium_noch_;
  bool p0_medium_; // also passes medium ID
  bool p0_tight_; // also passes tight ID
  unsigned p0_truth_; // 0 = default (assume jet), 1 = prompt photon, 2 = prompt electron
  bool p0_fsr_; // true if mother is a charged lepton, false otherwise

  // met
  float met_;
  float met_phi_;
  float tk_met_;
  float tk_met_phi_;
  float puppi_met_;
  float puppi_met_phi_;

  // composite variables
  float l0met_mt_;

  // composite variables defined inf n_p >= 1
  float l0p0_dr_;
  float l0p0_dphi_;
  float l0p0_M_;
  float reco_phi_;
  float reco_phi_f_;
  float reco_tk_phi_;
  float reco_tk_phi_f_;
  float reco_puppi_phi_;
  float reco_puppi_phi_f_;


  // vetos
  unsigned n_vm_; // number of additional veto muons
  float vm_p0_dr_; // if n_vm >= 1, deltaR between photon and closest veto muon

  // event weights
  float wt_def_; // default weight
  float wt_pf_; // prefiring weight
  float wt_pu_;
  float wt_l0_; // trk/ID/Iso weight for l0
  float wt_trg_l0_; // trigger weight for l0
  float wt_l1_; // trk/ID/Iso weight for m1
  float wt_p0_; // ID/iso weight for p0
  float wt_p0_fake_; // Photon fake factor
  float wt_p0_highpt_fake_; // Photon fake factor for high pT photons
  float wt_p0_e_fake_; // Electron -> photon fake factor

  // event weights for systematics (relative to nominal)
  float wt_sc_0_;
  float wt_sc_1_;
  float wt_sc_2_;
  float wt_sc_3_;
  float wt_sc_4_;
  float wt_sc_5_;

  std::vector<float> wt_pdf_;

  float wt_pu_hi_;
  float wt_pu_lo_;
  float wt_pf_hi_;
  float wt_pf_lo_;
  float wt_l0_hi_;
  float wt_l0_lo_;
  float wt_trg_l0_hi_;
  float wt_trg_l0_lo_;
  float wt_p0_hi_;
  float wt_p0_lo_;
  float wt_p0_fake_err_;
  int wt_p0_fake_bin_;
  float wt_p0_e_fake_hi_;
  float wt_p0_e_fake_lo_;

  // truth variables for Wgamma events
  bool is_wg_gen_;
  unsigned gen_nparts_;
  unsigned gen_pdgid_;
  bool gen_l0_match_;
  bool gen_p0_match_;
  float lhe_l0_pt_;
  float lhe_l0_eta_;
  float lhe_p0_pt_;
  float lhe_p0_eta_;
  float lhe_l0p0_dr_;
  float lhe_p0j_dr_;
  float lhe_j_pt_;
  bool lhe_frixione_;

  float gen_l0_pt_;
  float gen_l0_eta_;
  float gen_p0_pt_;
  float gen_p0_eta_;
  float gen_phi_;
  float gen_phi_f_;
  float true_phi_;
  float true_phi_f_;
  int gen_l0_q_;
  float gen_met_;
  float gen_met_phi_;
  float gen_l0p0_dr_;
  // float gen_dy_mll_;
  // float gen_n2_pt_;

  mutable TRandom3 rng;

  void PhotonIsoCorrector(ac::Photon *p, double rho);


 public:
  WGDataAnalysis(std::string const& name);
  virtual ~WGDataAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

  void SetDefaults();

  void CompressVars();
};
}  // namespace ac

#endif
