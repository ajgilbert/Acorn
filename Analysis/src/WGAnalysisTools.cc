#include "Acorn/Analysis/interface/WGAnalysisTools.h"
#include <algorithm>
#include <map>
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/NTupler/interface/city.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMath.h"
#include "boost/lexical_cast.hpp"

namespace ac {

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

}  // namespace ac
