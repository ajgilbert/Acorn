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
#include "Acorn/NTupler/interface/EventInfo.h"

namespace ac {

  WGAnalysis::WGAnalysis(std::string const& name)
      : ModuleBase(name), fs_(nullptr){}

  WGAnalysis::~WGAnalysis() { ; }

  int WGAnalysis::PreAnalysis() {
    if (fs_) {
      tree_ = fs_->make<TTree>("WGAnalysis", "WGAnalysis");
      tree_->Branch("phi1", &phi1_);
      tree_->Branch("rand_phi1", &rand_phi1_);
      tree_->Branch("w_pt", &w_pt_);
      tree_->Branch("g_pt", &g_pt_);
      tree_->Branch("n_pt", &n_pt_);
      tree_->Branch("l_pt", &l_pt_);
      tree_->Branch("l_eta", &l_eta_);
      tree_->Branch("type", &type_);
      tree_->Branch("nparts", &nparts_);
      tree_->Branch("weight_100004", &weight_100004_);
      // tree_->Branch("wt", &wt_);
      // tree_->Branch("n_jets", &n_jets_);
    }
    return 0;
  }

  int WGAnalysis::Execute(TreeEvent* event) {

    auto lhe_parts = event->GetPtrVec<ac::GenParticle>("lheParticles");
    nparts_ = lhe_parts.size();

    auto info = event->GetPtr<ac::EventInfo>("eventInfo");
    weight_100004_ = info->lheWeights().at(100001);

    ROOT::Math::PxPyPzEVector wg_system;
    ROOT::Math::PxPyPzEVector w_system;
    ROOT::Math::PxPyPzEVector charged_lepton;
    ROOT::Math::PxPyPzEVector neutrino;
    ROOT::Math::PxPyPzEVector photon;
    for (auto const& part : lhe_parts) {
      if (IsLepton(*part) || IsPhoton(*part)) {
        wg_system += part->vector();
        if (IsLepton(*part)) {
          w_system += part->vector();
          if (IsChargedLepton(*part)) {
            charged_lepton += part->vector();
          }
          if (IsNeutrino(*part)) {
            neutrino += part->vector();
          }
        }
        if (IsPhoton(*part)) {
          photon += part->vector();
        }
      }
    }
    w_pt_ = w_system.pt();
    g_pt_ = photon.pt();
    l_pt_ = charged_lepton.pt();
    l_eta_ = charged_lepton.eta();
    n_pt_ = neutrino.pt();
    // std::cout << "wg_system: " << wg_system << "\n";
    // std::cout << "boost: " << wg_system.BoostToCM() << "\n";
    auto boost = ROOT::Math::Boost(wg_system.BoostToCM());
    auto r_uvec = wg_system.Vect().unit();
    ac::GenParticle c_charged_lepton;
    ac::GenParticle c_neutrino;
    ac::GenParticle c_photon;

    for (auto const& part : lhe_parts) {
      if (IsLepton(*part) || IsPhoton(*part)) {
        auto boosted_part = ac::GenParticle(*part);
        boosted_part.setVector(boost(part->vector()));

        if (IsChargedLepton(*part)) {
          c_charged_lepton = boosted_part;
        }
        if (IsNeutrino(*part)) {
          c_neutrino = boosted_part;
        }
        if (IsPhoton(*part)) {
          c_photon = boosted_part;
          // boosted_part.Print();
        }

      }
    }
    ac::GenParticle c_w_boson = c_charged_lepton;
    c_w_boson.setVector(c_w_boson.vector() + c_neutrino.vector());
    auto z_uvec = c_w_boson.vector().Vect().unit();
    auto y_uvec = z_uvec.Cross(r_uvec).unit();
    auto x_uvec = y_uvec.Cross(z_uvec).unit();

    ROOT::Math::Rotation3D rotator(x_uvec, y_uvec, z_uvec);
    rotator.Invert();
    // std::cout << "c_w_boson: " << ROOT::Math::PxPyPzMVector(c_w_boson.vector()) << "\n";
    // std::cout << "c_photon:  " << ROOT::Math::PxPyPzMVector(c_photon.vector()) << "\n";

    ac::GenParticle r_w_boson = c_w_boson;
    r_w_boson.setVector(rotator(r_w_boson.vector()));

    ac::GenParticle r_charged_lepton = c_charged_lepton;
    r_charged_lepton.setVector(rotator(r_charged_lepton.vector()));

    ac::GenParticle r_neutrino = c_neutrino;
    r_neutrino.setVector(rotator(r_neutrino.vector()));

    ac::GenParticle r_photon = c_photon;
    r_photon.setVector(rotator(r_photon.vector()));

    // std::cout << "r_w_boson: " << ROOT::Math::PxPyPzMVector(r_w_boson.vector()) << "\n";
    // std::cout << "r_photon:  " << ROOT::Math::PxPyPzMVector(r_photon.vector()) << "\n";
    // std::cout << "r_rvec:    " << rotator(r_uvec) << "\n";

    // c_charged_lepton.Print();
    // c_neutrino.Print();
    if (r_charged_lepton.spin() * r_neutrino.spin() > 0) {
      type_ = 1;
      phi1_ = r_charged_lepton.phi();
      // std::cout << "Have W polarisation = 0\n";
    } else {
      type_ = 0;
      if (r_charged_lepton.spin() > 0) {
        if (r_charged_lepton.pdgId() > 0) type_ = 2;
        if (r_charged_lepton.pdgId() < 0) type_ = 3;
        phi1_ = r_charged_lepton.phi();
      }
      if (r_neutrino.spin() > 0) {
        phi1_ = r_neutrino.phi();
      }
    }


    if (rng.Uniform() > 0.5) {
      rand_phi1_ = ROOT::Math::VectorUtil::Phi_mpi_pi(ROOT::Math::Pi() - phi1_);
    }
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
