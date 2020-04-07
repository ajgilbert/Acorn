#ifndef Acorn_Analysis_WGAnalysisTools_h
#define Acorn_Analysis_WGAnalysisTools_h
#include <string>
#include <cstdint>
#include <algorithm>
#include "boost/range/algorithm_ext/erase.hpp"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include "TRandom3.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/Photon.h"

namespace ac {

struct WGSystem {
  ROOT::Math::PxPyPzEVector wg_system;
  ROOT::Math::PxPyPzEVector w_system;
  ROOT::Math::PxPyPzEVector charged_lepton;
  ROOT::Math::PxPyPzEVector true_neutrino;
  ROOT::Math::PxPyPzEVector neutrino;
  ROOT::Math::PxPyPzEVector photon;

  // Particles boosted into the CoM frame
  ac::Candidate c_charged_lepton;
  ac::Candidate c_neutrino;
  ac::Candidate c_photon;
  ac::Candidate c_w_boson;

  ac::Candidate r_w_boson;
  ac::Candidate r_charged_lepton;
  ac::Candidate r_neutrino;
  ac::Candidate r_photon;

  bool valid_mt;

  double Phi(unsigned lepton_charge);
  double SymPhi(unsigned lepton_charge);
};

WGSystem ProduceWGSystem(ac::Candidate const& lep, ac::Candidate const& neu,
                         ac::Candidate const& pho, bool reconstruct, TRandom3 & rng, bool verbose);

struct WGGenParticles {
  bool ok = true;
  unsigned nparts = 0;

  ac::GenParticle const* lhe_lep = nullptr;
  ac::GenParticle const* lhe_neu = nullptr;
  ac::GenParticle const* lhe_pho = nullptr;
  ac::GenParticle const* gen_lep = nullptr;
  ac::GenParticle const* gen_neu = nullptr;
  ac::GenParticle const* gen_pho = nullptr;

  std::vector<ac::GenParticle const*> viable_leptons;
  std::vector<ac::GenParticle const*> viable_photons;
};

WGGenParticles ProduceWGGenParticles(std::vector<GenParticle*> const& lhe_parts,
                                     std::vector<GenParticle*> const& gen_parts,
                                     double photon_dr = 0.7, unsigned version = 0);

bool IsElectron(ac::GenParticle const& p);
bool IsMuon(ac::GenParticle const& p);
bool IsTau(ac::GenParticle const& p);
bool IsChargedLepton(ac::GenParticle const& p);
bool IsNeutrino(ac::GenParticle const& p);
bool IsLepton(ac::GenParticle const& p);
bool IsPhoton(ac::GenParticle const& p);

std::vector<double> ExtractScaleVariations(ac::EventInfo const& info, int version);

}
#endif
