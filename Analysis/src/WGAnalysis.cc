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
      tree_->Branch("phi1", &phi1_);
      tree_->Branch("lhe_phi1", &lhe_phi1_);
      tree_->Branch("w_pt", &w_pt_);
      tree_->Branch("g_pt", &g_pt_);
      tree_->Branch("g_eta", &g_eta_);
      tree_->Branch("n_pt", &n_pt_);
      tree_->Branch("n_eta", &n_eta_);
      tree_->Branch("l_pt", &l_pt_);
      tree_->Branch("l_eta", &l_eta_);
      tree_->Branch("l_charge", &l_charge_);
      tree_->Branch("l_g_dr", &l_g_dr_);
      tree_->Branch("type", &type_);
      tree_->Branch("nparts", &nparts_);
      tree_->Branch("valid_mt", &valid_mt_);
      tree_->Branch("wt_C3w_0p0_", &wt_C3w_0p0_);
      tree_->Branch("wt_C3w_0p1_", &wt_C3w_0p1_);
      tree_->Branch("wt_C3w_0p2_", &wt_C3w_0p2_);
      tree_->Branch("wt_C3w_0p4_", &wt_C3w_0p4_);
      tree_->Branch("wt_C3w_1p0_", &wt_C3w_1p0_);
      // tree_->Branch("wt", &wt_);
      // tree_->Branch("n_jets", &n_jets_);
    }
    return 0;
  }

  int WGAnalysis::Execute(TreeEvent* event) {
    bool lhe_only = true;

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
    wt_C3w_1p0_ = 1.0;

    auto info = event->GetPtr<ac::EventInfo>("eventInfo");

    wt_C3w_0p0_ = info->lheWeights().at(100000);
    wt_C3w_0p1_ = info->lheWeights().at(100001);
    wt_C3w_0p2_ = info->lheWeights().at(100002);
    wt_C3w_0p4_ = info->lheWeights().at(100003);

    ac::GenParticle const* lhe_lep = nullptr;
    ac::GenParticle const* lhe_neu = nullptr;
    ac::GenParticle const* lhe_pho = nullptr;
    for (auto const& part : lhe_parts) {
      if (part->status() == 1) ++nparts_;
      if (IsChargedLepton(*part)) {
        lhe_lep = part;
      }
      if (IsNeutrino(*part)) {
        lhe_neu = part;
      }
      if (IsPhoton(*part)) {
        lhe_pho = part;
      }
    }

    if (std::abs(lhe_lep->pdgId()) == 15 && !lhe_only) {
      return 1;
    }

    ac::GenParticle const* gen_lep = nullptr;
    // ac::GenParticle const* gen_neu = nullptr;
    ac::GenParticle const* gen_pho = nullptr;

    std::vector<ac::GenParticle const*> viable_leptons;
    std::vector<ac::GenParticle const*> viable_photons;

    for (auto const& part : gen_parts) {
      // part->Print();
      if (IsChargedLepton(*part) && part->status() == 1) {
        viable_leptons.push_back(part);
        // gen_lep = part;
      }
      // if (IsNeutrino(*part)) {
      //   gen_neu = part;
      // }
      if (IsPhoton(*part) && part->status() == 1) {
        viable_photons.push_back(part);
      }
    }

    std::sort(viable_leptons.begin(), viable_leptons.end(), [&](ac::GenParticle const* c1, ac::GenParticle const* c2) {
      return ac::DeltaR(c1, lhe_lep) < ac::DeltaR(c2, lhe_lep);
    });
    std::sort(viable_photons.begin(), viable_photons.end(), [&](ac::GenParticle const* c1, ac::GenParticle const* c2) {
      return ac::DeltaR(c1, lhe_pho) < ac::DeltaR(c2, lhe_pho);
    });
    if (viable_leptons.size() && viable_photons.size()) {
      gen_lep = viable_leptons[0];
      gen_pho = viable_photons[0];
    } else if (!lhe_only) {
      std::cout << "No viable\n";
      return 1;
    }

    WGSystem lhe_sys = ProduceWGSystem(*lhe_lep, *lhe_neu, *lhe_pho, false, rng, false);
    WGSystem gen_sys;
    if (lhe_only) {
      gen_sys = ProduceWGSystem(*lhe_lep, *lhe_neu, *lhe_pho, true, rng, false);
    } else{
      gen_sys = ProduceWGSystem(*gen_lep, *gen_met, *gen_pho, true, rng, false);
    }

    w_pt_  = gen_sys.w_system.pt();
    g_pt_  = gen_sys.photon.pt();
    g_eta_ = gen_sys.photon.eta();
    l_pt_  = gen_sys.charged_lepton.pt();
    l_eta_ = gen_sys.charged_lepton.eta();
    l_charge_ = -1 * (lhe_lep->pdgId() / std::abs(lhe_lep->pdgId()));
    n_pt_  = gen_sys.neutrino.pt();
    n_eta_  = gen_sys.neutrino.eta();
    l_g_dr_ = ROOT::Math::VectorUtil::DeltaR(gen_sys.charged_lepton, gen_sys.photon);
    type_ = (lhe_lep->spin() * lhe_neu->spin()) > 0.;
    // std::cout << ">>>> Event\n";
    // lhe_lep->Print();
    // lhe_neu->Print();

    double lep_phi = gen_sys.r_charged_lepton.phi();
    // LHE particles don't have a charge assigned, use pdgid instead
    phi1_ = ROOT::Math::VectorUtil::Phi_mpi_pi(
        lhe_lep->pdgId() < 0 ? (lep_phi) : (lep_phi + ROOT::Math::Pi()));

    valid_mt_ = gen_sys.valid_mt;

    double lep_phi_lhe = lhe_sys.r_charged_lepton.phi();

    lhe_phi1_ = ROOT::Math::VectorUtil::Phi_mpi_pi(
        lhe_lep->pdgId() < 0 ? (lep_phi_lhe) : (lep_phi_lhe + ROOT::Math::Pi()));
    // std::cout << "Lepton charge is: " << l_charge_ << "\n";
    // std::cout << "phi_true = " << lhe_phi1_ << "\n";

    tree_->Fill();

    return 0;
  }
  int WGAnalysis::PostAnalysis() {
    return 0;
  }

  void WGAnalysis::PrintInfo() {}
}
