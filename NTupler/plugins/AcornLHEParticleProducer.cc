#include "Acorn/NTupler/plugins/AcornLHEParticleProducer.h"
#include <string>
#include <vector>
#include <bitset>
#include "TLorentzVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "Acorn/NTupler/interface/GenParticle.h"

// See: http://home.thep.lu.se/~leif/LHEF/classLHEF_1_1HEPEUP.html

AcornLHEParticleProducer::AcornLHEParticleProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::GenParticle>>(config),
      inputToken_(
          consumes<LHEEventProduct>(config.getParameter<edm::InputTag>("input"))),
      incomingP4Fix_(config.getParameter<bool>("incomingP4Fix")) {}

AcornLHEParticleProducer::~AcornLHEParticleProducer() { ; }

void AcornLHEParticleProducer::produce(edm::Event& event,
                                    const edm::EventSetup& setup) {
  edm::Handle<LHEEventProduct> lhe_handle;
  event.getByToken(inputToken_, lhe_handle);

  std::vector<lhef::HEPEUP::FiveVector> lhe_particles = lhe_handle->hepeup().PUP;

  output()->clear();
  output()->resize(lhe_particles.size(), ac::GenParticle());

  for (unsigned i = 0; i < lhe_handle->hepeup().PUP.size(); ++i) {
    ac::GenParticle& dest = output()->at(i);

    int status = lhe_handle->hepeup().ISTUP[i];

    ROOT::Math::PtEtaPhiMVector p4;

    if (incomingP4Fix_ && status == -1) {
      p4.SetPt(lhe_particles[i][2]);
      p4.SetM(lhe_particles[i][3]);
    } else {
      ROOT::Math::PxPyPzMVector p4c(lhe_particles[i][0], lhe_particles[i][1], lhe_particles[i][2],
                                   lhe_particles[i][4]);
      p4 = ROOT::Math::PtEtaPhiMVector(p4c);
    }

    dest.setVector(setVar("p4", p4));
    dest.setIndex(setVar("index", i));
    dest.setPdgId(setVar("pdgId", lhe_handle->hepeup().IDUP[i]));
    dest.setStatus(setVar("status", status));
    dest.setSpin(setVar("spin",lhe_handle->hepeup().SPINUP[i]));
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornLHEParticleProducer);
