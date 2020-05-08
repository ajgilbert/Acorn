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
#include "Acorn/NTupler/src/CMS_2020_PAS_SMP_20_005.h"

namespace ac {

  WGAnalysis::WGAnalysis(std::string const& name)
      : ModuleBase(name), fs_(nullptr), add_standalone_(false), add_rivet_(false) {}

  WGAnalysis::~WGAnalysis() { ; }

  int WGAnalysis::PreAnalysis() {
    if (fs_) {
      tree_ = fs_->make<TTree>("WGAnalysis", "WGAnalysis");

      tree_->Branch("is_wg_gen", &is_wg_gen_);
      tree_->Branch("gen_l0_pt", &gen_l0_pt_);
      tree_->Branch("gen_l0_eta", &gen_l0_eta_);
      tree_->Branch("gen_l0_q", &gen_l0_q_);
      tree_->Branch("gen_pdgid", &gen_pdgid_);
      tree_->Branch("gen_p0_pt", &gen_p0_pt_);
      tree_->Branch("gen_p0_eta", &gen_p0_eta_);
      tree_->Branch("gen_p0_frixione", &gen_p0_frixione_);
      tree_->Branch("gen_l0p0_dr", &gen_l0p0_dr_);
      tree_->Branch("gen_n0_pt", &gen_n0_pt_);
      tree_->Branch("gen_n0_eta", &gen_n0_eta_);
      tree_->Branch("gen_met_pt", &gen_met_pt_);
      tree_->Branch("gen_true_phi", &gen_true_phi_);
      tree_->Branch("gen_true_phi_f", &gen_true_phi_f_);


      if (add_rivet_) {
        tree_->Branch("is_wg_riv", &is_wg_riv_);
        tree_->Branch("riv_l0_pt", &riv_l0_pt_);
        tree_->Branch("riv_l0_eta", &riv_l0_eta_);
        tree_->Branch("riv_l0_q", &riv_l0_q_);
        tree_->Branch("riv_pdgid", &riv_pdgid_);
        tree_->Branch("riv_p0_pt", &riv_p0_pt_);
        tree_->Branch("riv_p0_eta", &riv_p0_eta_);
        tree_->Branch("riv_p0_frixione", &riv_p0_frixione_);
        tree_->Branch("riv_l0p0_dr", &riv_l0p0_dr_);
        tree_->Branch("riv_n0_pt", &riv_n0_pt_);
        tree_->Branch("riv_n0_eta", &riv_n0_eta_);
        tree_->Branch("riv_met_pt", &riv_met_pt_);
        tree_->Branch("riv_true_phi", &riv_true_phi_);
        tree_->Branch("riv_true_phi_f", &riv_true_phi_f_);
      }

      tree_->Branch("gen_phi", &gen_phi_);
      tree_->Branch("lhe_true_phi", &lhe_true_phi_);
      tree_->Branch("gen_phi_f", &gen_phi_f_);
      tree_->Branch("lhe_true_phi_f", &lhe_true_phi_f_);
      tree_->Branch("gen_W_pt", &gen_W_pt_);
      tree_->Branch("gen_l0p0_dphi", &gen_l0p0_dphi_);
      tree_->Branch("gen_Wp0_dphi", &gen_Wp0_dphi_);
      tree_->Branch("type", &type_);
      tree_->Branch("nparts", &nparts_);
      tree_->Branch("valid_mt", &valid_mt_);
      tree_->Branch("wt_def", &wt_def_);
      tree_->Branch("j1_pt", &j1_pt_);
      tree_->Branch("j1_eta", &j1_eta_);
      tree_->Branch("n_jets", &n_jets_);
      tree_->Branch("wt_C3w_0p0", &wt_C3w_0p0_);
      tree_->Branch("wt_C3w_0p1", &wt_C3w_0p1_);
      tree_->Branch("wt_C3w_0p2", &wt_C3w_0p2_);
      tree_->Branch("wt_C3w_0p4", &wt_C3w_0p4_);
      tree_->Branch("wt_C3w_0p67", &wt_C3w_0p67_);
      tree_->Branch("wt_C3w_1p0", &wt_C3w_1p0_);

      if (add_standalone_) {
        tree_->Branch("st_C3w_0p0", &st_C3w_0p0_);
        tree_->Branch("st_C3w_0p1", &st_C3w_0p1_);
        tree_->Branch("st_C3w_0p2", &st_C3w_0p2_);
        tree_->Branch("st_C3w_0p4", &st_C3w_0p4_);
        tree_->Branch("st_C3w_0p67", &st_C3w_0p67_);
        tree_->Branch("st_C3w_1p0", &st_C3w_1p0_);

        rw_ = std::make_shared<StandaloneReweight>("Acorn.Analysis.standalone_reweight", "rw_WG-EWDim6");
        replace_pdg_ = {
            {11, 13},
            {12, 14},
            {15, 13},
            {16, 14},
            {-11, -13},
            {-12, -14},
            {-15, -13},
            {-16, -14}
          };
      }
    }
    return 0;
  }

  int WGAnalysis::Execute(TreeEvent* event) {
    SetDefaults();
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

    // If this is the biased generation sample need to adjust the interpretation of the weights
    double rescale_st = 1.;
    if (info->lheWeights().count(100000)) {
      wt_C3w_0p0_ = info->lheWeights().at(100000) + 1.0;
      wt_C3w_0p1_ = info->lheWeights().at(100001) + 1.0;
      wt_C3w_0p2_ = info->lheWeights().at(100002) + 1.0;
      wt_C3w_0p4_ = info->lheWeights().at(100003) + 1.0;
      wt_C3w_0p67_ = info->lheWeights().at(100004) + 1.0;
      rescale_st = wt_C3w_0p0_;
    }

    if (add_standalone_) {
          std::vector<std::vector<double>> parts;
          std::vector<int> pdgs;
          std::vector<int> hels;
          std::vector<int> stats;

          for (auto const& p : lhe_parts) {
            if (std::abs(p->status()) != 1) {
              continue;
            }
            if (p->status() == -1) {
              parts.push_back({p->M(), 0., 0., p->pt()});
            }
            if (p->status() == 1) {
              parts.push_back({p->energy(), p->vector().px(), p->vector().py(), p->vector().pz()});
            }
            if (replace_pdg_.count(p->pdgId())) {
              pdgs.push_back(replace_pdg_[p->pdgId()]);
            } else {
              pdgs.push_back(p->pdgId());
            }
            hels.push_back(int(p->spin()));
            stats.push_back(p->status());
          }

          auto wts = rw_->ComputeWeights(parts, pdgs, hels, stats, info->lheAlphaS(), false, false);
          // std::cout << wts.size() << "\n";
          auto trans_wts = rw_->TransformWeights(wts);
          double param_card_Di = 0.1;
          // Not ideal, hardcoding 0.1 here
          double Ai = trans_wts[1] / param_card_Di;
          double Bii = trans_wts[2] / (param_card_Di * param_card_Di);
          st_C3w_0p0_ = 1.0;
          st_C3w_0p1_ = 1.0 + (1.5 * Ai) + (1.5 * 1.5 * Bii);
          st_C3w_0p2_ = 1.0 + (3.0 * Ai) + (3.0 * 3.0 * Bii);
          st_C3w_0p4_ = 1.0 + (6.0 * Ai) + (6.0 * 6.0 * Bii);
          st_C3w_0p67_ = 1.0 + (10.0 * Ai) + (10.0 * 10.0 * Bii);
          st_C3w_1p0_ = 1.0 + (15.0 * Ai) + (15.0 * 15.0 * Bii);
          st_C3w_0p0_ *= rescale_st;
          st_C3w_0p1_ *= rescale_st;
          st_C3w_0p2_ *= rescale_st;
          st_C3w_0p4_ *= rescale_st;
          st_C3w_0p67_ *= rescale_st;
          st_C3w_1p0_ *= rescale_st;

          // std::cout << "0p0: " << wt_C3w_0p0_ << "\t" << st_C3w_0p0_ << "\t"  << (st_C3w_0p0_ / wt_C3w_0p0_) << "\n";
          // std::cout << "0p1: " << wt_C3w_0p1_ << "\t" << st_C3w_0p1_ << "\t"  << (st_C3w_0p1_ / wt_C3w_0p1_) << "\n";
          // std::cout << "0p2: " << wt_C3w_0p2_ << "\t" << st_C3w_0p2_ << "\t"  << (st_C3w_0p2_ / wt_C3w_0p2_) << "\n";
          // std::cout << "0p4: " << wt_C3w_0p4_ << "\t" << st_C3w_0p4_ << "\t"  << (st_C3w_0p4_ / wt_C3w_0p4_) << "\n";
          // std::cout << "0p67: " << wt_C3w_0p67_ << "\t" << st_C3w_0p67_ << "\t"  << (st_C3w_0p67_ / wt_C3w_0p67_) << "\n";
    }

    WGGenParticles parts = ProduceWGGenParticles(lhe_parts, gen_parts);
    // WGGenParticles parts = ProduceWGGenParticles(lhe_parts, gen_parts, 0.7, 1);
    nparts_ = parts.nparts;

    is_wg_gen_ = parts.ok;

    if (is_wg_gen_) {
      WGSystem lhe_sys = ProduceWGSystem(*parts.lhe_lep, *parts.lhe_neu, *parts.lhe_pho, false, rng, false);
      WGSystem gen_sys;
      if (lhe_only) {
        gen_sys = ProduceWGSystem(*parts.lhe_lep, *parts.lhe_neu, *parts.lhe_pho, true, rng, false);
      } else{
        gen_sys = ProduceWGSystem(*parts.gen_lep, *gen_met, *parts.gen_pho, true, rng, false);
      }
      WGSystem gen_true_sys = ProduceWGSystem(*parts.gen_lep, *parts.gen_neu, *parts.gen_pho, false, rng, false);

      gen_W_pt_  = gen_sys.w_system.pt();
      gen_p0_pt_  = gen_sys.photon.pt();
      gen_p0_eta_ = gen_sys.photon.eta();
      gen_p0_frixione_ = ac::FrixioneIso(*parts.lhe_pho, lhe_parts, 0.4);
      gen_l0_pt_  = gen_sys.charged_lepton.pt();
      gen_l0_eta_ = gen_sys.charged_lepton.eta();
      gen_pdgid_ = std::abs(parts.lhe_lep->pdgId());
      gen_l0_q_ = -1 * (parts.lhe_lep->pdgId() / std::abs(parts.lhe_lep->pdgId()));
      gen_n0_pt_  = gen_sys.neutrino.pt();
      gen_n0_eta_  = gen_sys.neutrino.eta();
      gen_met_pt_ = gen_met->pt();
      gen_l0p0_dr_ = ROOT::Math::VectorUtil::DeltaR(gen_sys.charged_lepton, gen_sys.photon);
      gen_l0p0_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(gen_sys.charged_lepton, gen_sys.photon);
      gen_Wp0_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(gen_sys.w_system, gen_sys.photon);
      type_ = (parts.lhe_lep->spin() * parts.lhe_neu->spin()) > 0.;

      // Define the gen_phi (using final state GenParticles + MET)
      gen_phi_ = gen_sys.Phi(parts.gen_lep->pdgId() < 0);
      gen_phi_f_ = gen_sys.SymPhi(parts.gen_lep->pdgId() < 0);

      // Define the the true_phi (using final state GenParticles)
      gen_true_phi_ = gen_true_sys.Phi(parts.gen_lep->pdgId() < 0);
      gen_true_phi_f_ = gen_true_sys.SymPhi(parts.gen_lep->pdgId() < 0);

      valid_mt_ = gen_sys.valid_mt;

      // Define the lhe_true_phi (using LHE particles)
      lhe_true_phi_ = lhe_sys.Phi(parts.lhe_lep->pdgId() < 0);
      lhe_true_phi_f_ = lhe_sys.SymPhi(parts.lhe_lep->pdgId() < 0);
    }

    if (add_rivet_) {
      WGammaRivetVariables const* rivet = event->GetPtr<WGammaRivetVariables>("rivetVariables");
      is_wg_riv_ = rivet->is_wg_gen;
      riv_l0_pt_ = rivet->l0_pt;
      riv_l0_eta_ =rivet->l0_eta;
      riv_l0_q_ = rivet->l0_q;
      riv_pdgid_ = rivet->l0_abs_pdgid;
      riv_p0_pt_ = rivet->p0_pt;
      riv_p0_eta_ = rivet->p0_eta;
      riv_p0_frixione_ = rivet->p0_frixione;
      riv_l0p0_dr_ = rivet->l0p0_dr;
      riv_n0_pt_ = rivet->n0_pt;
      riv_n0_eta_ = rivet->n0_eta;
      riv_met_pt_ = rivet->met_pt;
      riv_true_phi_ = rivet->true_phi;
      riv_true_phi_f_ = rivet->true_phi_f;
    }


    tree_->Fill();

    return 0;
  }
  int WGAnalysis::PostAnalysis() {
    return 0;
  }

  void WGAnalysis::PrintInfo() {}


  void WGAnalysis::SetDefaults() {
    is_wg_gen_ = false;
    gen_l0_pt_ = 0.;
    gen_l0_eta_ = 0.;
    gen_l0_q_ = 0;
    gen_pdgid_ = 0;
    gen_p0_pt_ = 0.;
    gen_p0_eta_ = 0.;
    gen_p0_frixione_ = false;
    gen_l0p0_dr_ = 0.;
    gen_n0_pt_ = 0.;
    gen_n0_eta_ = 0.;
    gen_met_pt_ = 0.;
    gen_true_phi_ = 0.;
    gen_true_phi_f_ = 0.;

    is_wg_riv_ = false;
    riv_l0_pt_ = 0.;
    riv_l0_eta_ = 0.;
    riv_l0_q_ = 0;
    riv_pdgid_ = 0;
    riv_p0_pt_ = 0.;
    riv_p0_eta_ = 0.;
    riv_p0_frixione_ = false;
    riv_l0p0_dr_ = 0.;
    riv_n0_pt_ = 0.;
    riv_n0_eta_ = 0.;
    riv_met_pt_ = 0.;
    riv_true_phi_ = 0.;
    riv_true_phi_f_ = 0.;
    gen_phi_ = 0.;
    gen_phi_f_ = 0.;
    lhe_true_phi_ = 0.;
    lhe_true_phi_f_ = 0.;
    type_ = 0;
    wt_def_ = 1.;
    wt_C3w_0p0_ = 1.;
    wt_C3w_0p1_ = 1.;
    wt_C3w_0p2_ = 1.;
    wt_C3w_0p4_ = 1.;
    wt_C3w_0p67_ = 1.;
    wt_C3w_1p0_ = 1.;
    gen_W_pt_ = 0.;
    gen_l0p0_dphi_ = 0.;
    gen_Wp0_dphi_ = 0.;
    j1_pt_ = 0.;
    j1_eta_ = 0.;
    n_jets_ = 0;
    nparts_ = 0;
    valid_mt_ = false;

    // Weights computed using standalone ME RW
    st_C3w_0p0_ = 1.;
    st_C3w_0p1_ = 1.;
    st_C3w_0p2_ = 1.;
    st_C3w_0p4_ = 1.;
    st_C3w_0p67_ = 1.;
    st_C3w_1p0_ = 1.;
  }

}

