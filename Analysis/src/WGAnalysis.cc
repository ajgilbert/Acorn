#include <algorithm>
#include <map>
#include "TMath.h"
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"
#include "RooRealVar.h"
#include "Acorn/Analysis/interface/WGAnalysis.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/Analysis/interface/WGAnalysisTools.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

namespace ac {

  WGAnalysis::WGAnalysis(std::string const& name)
      : ModuleBase(name), fs_(nullptr){}

  WGAnalysis::~WGAnalysis() { ; }

  int WGAnalysis::PreAnalysis() {
    if (fs_) {
      tree_ = fs_->make<TTree>("WGAnalysis", "WGAnalysis");
      tree_->Branch("gen_phi", &gen_phi_);
      tree_->Branch("true_phi", &true_phi_);
      tree_->Branch("lhe_true_phi", &lhe_true_phi_);
      tree_->Branch("gen_phi_f", &gen_phi_f_);
      tree_->Branch("true_phi_f", &true_phi_f_);
      tree_->Branch("lhe_true_phi_f", &lhe_true_phi_f_);
      tree_->Branch("w_pt", &w_pt_);
      tree_->Branch("g_pt", &g_pt_);
      tree_->Branch("g_eta", &g_eta_);
      tree_->Branch("n_pt", &n_pt_);
      tree_->Branch("n_eta", &n_eta_);
      tree_->Branch("l_pt", &l_pt_);
      tree_->Branch("l_eta", &l_eta_);
      tree_->Branch("l_charge", &l_charge_);
      tree_->Branch("l_g_dr", &l_g_dr_);
      tree_->Branch("l_g_dphi", &l_g_dphi_);
      tree_->Branch("w_g_dphi", &w_g_dphi_);
      tree_->Branch("type", &type_);
      tree_->Branch("nparts", &nparts_);
      tree_->Branch("valid_mt", &valid_mt_);
      tree_->Branch("wt_def", &wt_def_);
      tree_->Branch("wt_C3w_0p0", &wt_C3w_0p0_);
      tree_->Branch("wt_C3w_0p1", &wt_C3w_0p1_);
      tree_->Branch("wt_C3w_0p2", &wt_C3w_0p2_);
      tree_->Branch("wt_C3w_0p4", &wt_C3w_0p4_);
      tree_->Branch("wt_C3w_0p67", &wt_C3w_0p67_);
      tree_->Branch("wt_C3w_1p0", &wt_C3w_1p0_);
    }
    return 0;
  }

  int WGAnalysis::Execute(TreeEvent* event) {
    bool lhe_only = false;

    auto lhe_parts = event->GetPtrVec<ac::GenParticle>("lheParticles");
    std::vector<GenParticle *> gen_parts;
    Met* gen_met;
    if (!lhe_only) {
      gen_parts = event->GetPtrVec<ac::GenParticle>("genParticles");
      gen_met = event->GetPtrVec<ac::Met>("genMet")[0];
    }
    nparts_ = 0;

    wt_C3w_0p0_ = 1.0;
    wt_C3w_0p1_ = 1.0;
    wt_C3w_0p2_ = 1.0;
    wt_C3w_0p4_ = 1.0;
    wt_C3w_0p67_ = 1.0;
    wt_C3w_1p0_ = 1.0;

    auto info = event->GetPtr<ac::EventInfo>("eventInfo");

    wt_def_ = info->totalWeight();

    wt_C3w_0p0_ = info->lheWeights().at(100000);
    wt_C3w_0p1_ = info->lheWeights().at(100001);
    wt_C3w_0p2_ = info->lheWeights().at(100002);
    wt_C3w_0p4_ = info->lheWeights().at(100003);
    wt_C3w_0p67_ = info->lheWeights().at(100004);

    WGGenParticles parts = ProduceWGGenParticles(lhe_parts, gen_parts);
    nparts_ = parts.nparts;

    if (!parts.ok) {
      return 1;
    }

    WGSystem lhe_sys = ProduceWGSystem(*parts.lhe_lep, *parts.lhe_neu, *parts.lhe_pho, false, rng, false);
    WGSystem gen_sys;
    if (lhe_only) {
      gen_sys = ProduceWGSystem(*parts.lhe_lep, *parts.lhe_neu, *parts.lhe_pho, true, rng, false);
    } else{
      gen_sys = ProduceWGSystem(*parts.gen_lep, *gen_met, *parts.gen_pho, true, rng, false);
    }
    WGSystem gen_true_sys = ProduceWGSystem(*parts.gen_lep, *parts.gen_neu, *parts.gen_pho, false, rng, false);

    w_pt_  = gen_sys.w_system.pt();
    g_pt_  = gen_sys.photon.pt();
    g_eta_ = gen_sys.photon.eta();
    l_pt_  = gen_sys.charged_lepton.pt();
    l_eta_ = gen_sys.charged_lepton.eta();
    l_charge_ = -1 * (parts.lhe_lep->pdgId() / std::abs(parts.lhe_lep->pdgId()));
    n_pt_  = gen_sys.neutrino.pt();
    n_eta_  = gen_sys.neutrino.eta();
    l_g_dr_ = ROOT::Math::VectorUtil::DeltaR(gen_sys.charged_lepton, gen_sys.photon);
    l_g_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(gen_sys.charged_lepton, gen_sys.photon);
    w_g_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(gen_sys.w_system, gen_sys.photon);
    type_ = (parts.lhe_lep->spin() * parts.lhe_neu->spin()) > 0.;

    // Define the gen_phi (using final state GenParticles + MET)
    gen_phi_ = gen_sys.Phi(parts.gen_lep->pdgId() < 0);
    gen_phi_f_ = gen_sys.SymPhi(parts.gen_lep->pdgId() < 0);

    // Define the the true_phi (using final state GenParticles)
    true_phi_ = gen_true_sys.Phi(parts.gen_lep->pdgId() < 0);
    true_phi_f_ = gen_true_sys.SymPhi(parts.gen_lep->pdgId() < 0);

    valid_mt_ = gen_sys.valid_mt;

    // Define the lhe_true_phi (using LHE particles)
    lhe_true_phi_ = lhe_sys.Phi(parts.lhe_lep->pdgId() < 0);
    lhe_true_phi_f_ = lhe_sys.SymPhi(parts.lhe_lep->pdgId() < 0);

    tree_->Fill();

    return 0;
  }
  int WGAnalysis::PostAnalysis() {
    return 0;
  }

  void WGAnalysis::PrintInfo() {}
}
