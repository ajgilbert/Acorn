#include "Acorn/NTupler/plugins/AcornGenParticleProducer.h"
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

AcornGenParticleProducer::AcornGenParticleProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::GenParticle>>(config),
      inputToken_(
          consumes<edm::View<reco::GenParticle>>(config.getParameter<edm::InputTag>("input"))) {}

AcornGenParticleProducer::~AcornGenParticleProducer() { ; }

void AcornGenParticleProducer::produce(edm::Event& event,
                                    const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::GenParticle> > parts_handle;
  event.getByToken(inputToken_, parts_handle);

  output()->clear();
  output()->resize(parts_handle->size(), ac::GenParticle());

  for (unsigned i = 0; i < parts_handle->size(); ++i) {
    reco::GenParticle const& src = parts_handle->at(i);
    ac::GenParticle& dest = output()->at(i);

    dest.setVector(setVar("p4", src.polarP4()));
    dest.setCharge(setVar("charge", src.charge()));
    dest.setIndex(setVar("index", static_cast<int>(parts_handle->refAt(i).key())));
    dest.setPdgId(setVar("pdgId", src.pdgId()));
    dest.setStatus(setVar("status", src.status()));
    dest.setStatusFlags(setVar("statusFlags", src.statusFlags()));

    std::vector<int> mothers(src.motherRefVector().size(), 0);
    for (unsigned j = 0; j < src.motherRefVector().size(); ++j) {
      mothers[j] = static_cast<int>(src.motherRefVector().at(j).key());
    }
    dest.setMothers(setVar("mothers", mothers));

    std::vector<int> daughters(src.daughterRefVector().size(), 0);
    for (unsigned j = 0; j < src.daughterRefVector().size(); ++j) {
      daughters[j] = static_cast<int>(src.daughterRefVector().at(j).key());
    }
    dest.setDaughters(setVar("daughters", daughters));
    // dest.Print();
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornGenParticleProducer);
