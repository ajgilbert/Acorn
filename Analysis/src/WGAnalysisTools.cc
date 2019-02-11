#include "Acorn/Analysis/interface/WGAnalysisTools.h"
#include <algorithm>
#include <map>
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/city.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMath.h"
#include "boost/lexical_cast.hpp"

namespace ac {

double WGSystem::Phi(unsigned lepton_charge) {
  double lep_phi = r_charged_lepton.phi();
  return ROOT::Math::VectorUtil::Phi_mpi_pi(lepton_charge > 0 ? (lep_phi)
                                                              : (lep_phi + ROOT::Math::Pi()));
}
double WGSystem::SymPhi(unsigned lepton_charge) {
  double lep_phi = Phi(lepton_charge);
  if (lep_phi > ROOT::Math::Pi() / 2.) {
    return ROOT::Math::Pi() - lep_phi;
  } else if (lep_phi < -1. * (ROOT::Math::Pi() / 2.)) {
    return -1. * (ROOT::Math::Pi() + lep_phi);
  } else {
    return lep_phi;
  }
}

  WGSystem ProduceWGSystem(ac::Candidate const& lep, ac::Candidate const& neu,
                                       ac::Candidate const& pho, bool reconstruct, TRandom3 & rng, bool verbose) {
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

  WGGenParticles ProduceWGGenParticles(std::vector<GenParticle*> const& lhe_parts,
                                       std::vector<GenParticle*> const& gen_parts) {
    WGGenParticles info;

    for (auto const& part : lhe_parts) {
      if (part->status() == 1) ++info.nparts;
      if (IsChargedLepton(*part)) {
        info.lhe_lep = part;
      }
      if (IsNeutrino(*part)) {
        info.lhe_neu = part;
      }
      if (IsPhoton(*part)) {
        info.lhe_pho = part;
      }
    }

    if (std::abs(info.lhe_lep->pdgId()) == 15) {
      info.ok = false;
      return info;
    }

    for (auto const& part : gen_parts) {
      if (IsChargedLepton(*part) && part->statusFlags().isPrompt() && part->status() == 1) {
        info.viable_leptons.push_back(part);
        // gen_lep = part;
      }
      if (IsNeutrino(*part) && part->statusFlags().isPrompt() && part->status() == 1) {
        if (info.gen_neu) {
          std::cout << "Found more than one neutrino!\n";
          info.gen_neu->Print();
          part->Print();
        }
        info.gen_neu = part;
      }
      // Use an arbitrary pT > 1 GeV cut to throw away some FSR noise
      if (IsPhoton(*part) && part->statusFlags().isPrompt() && part->status() == 1 && part->pt() > 1.0) {
        info.viable_photons.push_back(part);
      }
    }

    std::sort(info.viable_leptons.begin(), info.viable_leptons.end(), [&](ac::GenParticle const* c1, ac::GenParticle const* c2) {
      return ac::DeltaR(c1, info.lhe_lep) < ac::DeltaR(c2, info.lhe_lep);
    });
    std::sort(info.viable_photons.begin(), info.viable_photons.end(), [&](ac::GenParticle const* c1, ac::GenParticle const* c2) {
      return c1->pt() > c2->pt();
    });
    if (info.viable_leptons.size() && info.viable_photons.size() && info.gen_neu) {
      info.gen_lep = info.viable_leptons[0];
      info.gen_pho = info.viable_photons[0];
    } else {
      // std::cout << "Not viable:" << info.viable_leptons.size() << "\t" << info.viable_photons.size() << "\t" << info.gen_neu << "\n";
      info.ok = false;
    }

    // if (info.viable_photons.size() == 0) {
    //     std::cout << "LHE particles:\n";
    //     for (auto const* p : lhe_parts) p->Print();
    //     std::cout << "GEN particles:\n";
    //     for (auto const* p : gen_parts) p->Print();
    // }

    return info;
  }

  bool IsElectron(ac::GenParticle const& p) {
    return  std::abs(p.pdgId()) == 11;
  }

  bool IsMuon(ac::GenParticle const& p) {
    return  std::abs(p.pdgId()) == 13;
  }
  bool IsTau(ac::GenParticle const& p) {
    return  std::abs(p.pdgId()) == 15;
  }

  bool IsChargedLepton(ac::GenParticle const& p) {
    return std::abs(p.pdgId()) == 11 || std::abs(p.pdgId()) == 13 || std::abs(p.pdgId()) == 15;
  }
  bool IsNeutrino(ac::GenParticle const& p) {
    return std::abs(p.pdgId()) == 12 || std::abs(p.pdgId()) == 14 || std::abs(p.pdgId()) == 16;
  }
  bool IsLepton(ac::GenParticle const& p) {
    return IsChargedLepton(p) || IsNeutrino(p);
  }

  bool IsPhoton(ac::GenParticle const& p) {
    return p.pdgId() == 22;
  }

}  // namespace ac
