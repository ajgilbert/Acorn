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
      tree_->Branch("weight_C3w_0p1", &weight_C3w_0p1_);
      tree_->Branch("weight_C3w_0p2", &weight_C3w_0p2_);
      tree_->Branch("weight_C3w_0p4", &weight_C3w_0p4_);
      tree_->Branch("weight_C3w_0p8", &weight_C3w_0p8_);
      tree_->Branch("weight_C3w_1p6", &weight_C3w_1p6_);
      // tree_->Branch("wt", &wt_);
      // tree_->Branch("n_jets", &n_jets_);
    }
    return 0;
  }

  WGSystem WGAnalysis::ProduceWGSystem(ac::Candidate const& lep, ac::Candidate const& neu,
                                       ac::Candidate const& pho, bool reconstruct, bool verbose) const {
    WGSystem sys;

    sys.charged_lepton += lep.vector();
    sys.true_neutrino += neu.vector();
    sys.photon += pho.vector();


    if (reconstruct) {
      // Interpret neu as the MET
      ROOT::Math::PxPyPzEVector met(
          ROOT::Math::PxPyPzMVector(neu.vector().px(), neu.vector().py(), 0., 0.));
      double mt = std::sqrt(2. * lep.pt() * met.pt() * (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(lep.vector(), met))));
      double w_mass = 80.385;

      // std::cout << mt << "\t" << w_mass << "\n";

      // Implement eqn 10 of CERN-TH-2017-185
      double neu_eta = 0.;
      if (mt < w_mass) {
        double delta2 = (w_mass * w_mass - mt * mt) / (2. * lep.pt() * met.pt());
        double delta = std::sqrt(delta2);
        double delta_eta = std::log(1. + (delta * std::sqrt(2. + delta2)) + delta2);
        if (rng.Uniform() > 0.5) {
          delta_eta *= -1.;
        }
        neu_eta = lep.eta() + delta_eta;
        // double neu_eta_alt = lep.eta() - delta_eta;
        // ROOT::Math::PxPyPzEVector rec_neu_alt(ROOT::Math::PtEtaPhiMVector(met.pt(), neu_eta_alt, met.phi(), 0.));
        // std::cout << "True neutrino: " << ROOT::Math::PxPyPzEVector(neu.vector()) << "\n";
        // std::cout << "Reco neutrino: " << rec_neu << "\n";
        // std::cout << "Alt neutrino:  " << rec_neu_alt << "\n";
        sys.valid_mt = true;
      } else {
        neu_eta = lep.eta();
        sys.valid_mt = false;
        // std::cout << "mT is greater than mW!\n";
        // std::cout << "True neutrino: " << ROOT::Math::PxPyPzEVector(neu.vector()) << "\n";
        // std::cout << "Reco neutrino: " << rec_neu << "\n";
      }
      ROOT::Math::PxPyPzEVector rec_neu(ROOT::Math::PtEtaPhiMVector(met.pt(), neu_eta, met.phi(), 0.));
      sys.neutrino = rec_neu;
    } else {
      sys.neutrino += neu.vector();
    }

    sys.wg_system += sys.charged_lepton;
    sys.wg_system += sys.neutrino;
    sys.wg_system += sys.photon;

    sys.w_system += sys.charged_lepton;
    sys.w_system += sys.neutrino;



    auto boost = ROOT::Math::Boost(sys.wg_system.BoostToCM());

    auto r_uvec = sys.wg_system.Vect().unit();

    sys.c_charged_lepton = lep;
    sys.c_neutrino = neu;
    sys.c_neutrino.setVector(ROOT::Math::PtEtaPhiMVector(sys.neutrino));
    sys.c_photon = pho;

    sys.c_charged_lepton.setVector(boost(sys.c_charged_lepton.vector()));
    sys.c_neutrino.setVector(boost(sys.c_neutrino.vector()));
    sys.c_photon.setVector(boost(sys.c_photon.vector()));

    sys.c_w_boson = sys.c_charged_lepton;
    sys.c_w_boson.setVector(sys.c_w_boson.vector() + sys.c_neutrino.vector());

    auto z_uvec = sys.c_w_boson.vector().Vect().unit();
    auto y_uvec = z_uvec.Cross(r_uvec).unit();
    auto x_uvec = y_uvec.Cross(z_uvec).unit();

    ROOT::Math::Rotation3D rotator(x_uvec, y_uvec, z_uvec);
    rotator.Invert();
    // std::cout << "c_w_boson: " << ROOT::Math::PxPyPzMVector(c_w_boson.vector()) << "\n";
    // std::cout << "c_photon:  " << ROOT::Math::PxPyPzMVector(c_photon.vector()) << "\n";

    sys.r_w_boson = sys.c_w_boson;
    sys.r_w_boson.setVector(rotator(sys.r_w_boson.vector()));

    sys.r_charged_lepton = sys.c_charged_lepton;
    sys.r_charged_lepton.setVector(rotator(sys.r_charged_lepton.vector()));

    sys.r_neutrino = sys.c_neutrino;
    sys.r_neutrino.setVector(rotator(sys.r_neutrino.vector()));

    sys.r_photon = sys.c_photon;
    sys.r_photon.setVector(rotator(sys.r_photon.vector()));

    if (verbose) {
      std::cout << "================================\n";
      std::cout << "Lab frame (pt, eta, phi, M) / (px, py, pz, E):\n";
      std::cout << "Electron:" << ROOT::Math::PtEtaPhiMVector(sys.charged_lepton) << " / " << ROOT::Math::PxPyPzEVector(sys.charged_lepton) << "\n";
      std::cout << "Neutrino:" << ROOT::Math::PtEtaPhiMVector(sys.neutrino) << " / " << ROOT::Math::PxPyPzEVector(sys.neutrino) << "\n";
      std::cout << "W boson: " << ROOT::Math::PtEtaPhiMVector(sys.w_system) << " / " << ROOT::Math::PxPyPzEVector(sys.w_system) << "\n";
      std::cout << "Photon:  " << ROOT::Math::PtEtaPhiMVector(sys.photon) << " / " << ROOT::Math::PxPyPzEVector(sys.photon) << "\n";
      std::cout << "---------------------------------\n";
      std::cout << "Boosted to WG system (pt, eta, phi, M) / (px, py, pz, E):\n";
      std::cout << "Electron:" << ROOT::Math::PtEtaPhiMVector(sys.c_charged_lepton.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.c_charged_lepton.vector())  << "\n";
      std::cout << "Neutrino:" << ROOT::Math::PtEtaPhiMVector(sys.c_neutrino.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.c_neutrino.vector())  << "\n";
      std::cout << "W boson: " << ROOT::Math::PtEtaPhiMVector(sys.c_w_boson.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.c_w_boson.vector()) << "\n";
      std::cout << "Photon:  " << ROOT::Math::PtEtaPhiMVector(sys.c_photon.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.c_photon.vector())  << "\n";
      std::cout << "---------------------------------\n";
      std::cout << "After rotation (pt, eta, phi, M) / (px, py, pz, E):\n";
      std::cout << "Electron:" << ROOT::Math::PtEtaPhiMVector(sys.r_charged_lepton.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.r_charged_lepton.vector())  << "\n";
      std::cout << "Neutrino:" << ROOT::Math::PtEtaPhiMVector(sys.r_neutrino.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.r_neutrino.vector())  << "\n";
      std::cout << "W boson: " << ROOT::Math::PtEtaPhiMVector(sys.r_w_boson.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.r_w_boson.vector()) << "\n";
      std::cout << "Photon:  " << ROOT::Math::PtEtaPhiMVector(sys.r_photon.vector()) << " / " << ROOT::Math::PxPyPzEVector(sys.r_photon.vector())  << "\n";
      std::cout << "---------------------------------\n";
    }
    return sys;
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

    weight_C3w_0p1_ = 1.0;
    weight_C3w_0p2_ = 1.0;
    weight_C3w_0p4_ = 1.0;
    weight_C3w_0p8_ = 1.0;
    weight_C3w_1p6_ = 1.0;

    auto info = event->GetPtr<ac::EventInfo>("eventInfo");
    if (info->lheWeights().count(100000)) {
      weight_C3w_0p1_ = info->lheWeights().at(100000);
      weight_C3w_0p4_ = info->lheWeights().at(100002);
      weight_C3w_0p8_ = info->lheWeights().at(100003);
      weight_C3w_1p6_ = info->lheWeights().at(100004);
    };
    if (info->lheWeights().count(100001)) {
      weight_C3w_0p2_ = info->lheWeights().at(100001);
    }

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

    WGSystem lhe_sys = ProduceWGSystem(*lhe_lep, *lhe_neu, *lhe_pho, false, false);
    WGSystem gen_sys;
    if (lhe_only) {
      gen_sys = ProduceWGSystem(*lhe_lep, *lhe_neu, *lhe_pho, true, false);
    } else{
      gen_sys = ProduceWGSystem(*gen_lep, *gen_met, *gen_pho, true, false);
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

  bool WGAnalysis::IsChargedLepton(ac::GenParticle const& p) const {
    return std::abs(p.pdgId()) == 11 || std::abs(p.pdgId()) == 13 || std::abs(p.pdgId()) == 15;
  }
  bool WGAnalysis::IsNeutrino(ac::GenParticle const& p) const {
    return std::abs(p.pdgId()) == 12 || std::abs(p.pdgId()) == 14 || std::abs(p.pdgId()) == 16;
  }
  bool WGAnalysis::IsLepton(ac::GenParticle const& p) const {
    return IsChargedLepton(p) || IsNeutrino(p);
  }

  bool WGAnalysis::IsPhoton(ac::GenParticle const& p) const {
    return p.pdgId() == 22;
  }

}
