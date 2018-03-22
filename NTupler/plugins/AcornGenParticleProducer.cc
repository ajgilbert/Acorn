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
      input_(config.getParameter<edm::InputTag>("input")){
  consumes<edm::View<reco::GenParticle>>(input_);
}

AcornGenParticleProducer::~AcornGenParticleProducer() { ; }

void AcornGenParticleProducer::produce(edm::Event& event,
                                    const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::GenParticle> > parts_handle;
  event.getByLabel(input_, parts_handle);

  output()->clear();
  output()->resize(parts_handle->size(), ac::GenParticle());

  for (unsigned i = 0; i < parts_handle->size(); ++i) {
    reco::GenParticle const& src = parts_handle->at(i);
    ac::GenParticle& dest = output()->at(i);

    dest.set_vector(setVar("p4", src.polarP4()));
    dest.set_charge(setVar("charge", src.charge()));
    dest.set_index(setVar("index", static_cast<int>(parts_handle->refAt(i).key())));
    dest.set_pdgid(setVar("pdgId", src.pdgId()));
    dest.set_status(setVar("status", src.status()));
    dest.set_statusFlags(setVar("statusFlags", src.statusFlags()));

    std::vector<int> mothers(src.motherRefVector().size(), 0);
    for (unsigned j = 0; j < src.motherRefVector().size(); ++j) {
      mothers[j] = static_cast<int>(src.motherRefVector().at(j).key());
    }
    dest.set_mothers(setVar("mothers", mothers));

    std::vector<int> daughters(src.daughterRefVector().size(), 0);
    for (unsigned j = 0; j < src.daughterRefVector().size(); ++j) {
      daughters[j] = static_cast<int>(src.daughterRefVector().at(j).key());
    }
    dest.set_daughters(setVar("daughters", daughters));
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornGenParticleProducer);
