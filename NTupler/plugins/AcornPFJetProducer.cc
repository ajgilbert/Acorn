#include "Acorn/NTupler/plugins/AcornPFJetProducer.h"
#include <string>
#include <vector>
#include "Acorn/NTupler/interface/PFJet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

AcornPFJetProducer::AcornPFJetProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::PFJet>>(config),
      inputToken_(consumes<edm::View<pat::Jet>>(config.getParameter<edm::InputTag>("input"))),
      jetidSelector_(config.getParameter<edm::ParameterSet>("jetID")) {}

AcornPFJetProducer::~AcornPFJetProducer() { ; }

void AcornPFJetProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<pat::Jet>> jets_handle;
  event.getByToken(inputToken_, jets_handle);

  output()->clear();
  output()->resize(jets_handle->size(), ac::PFJet());

  for (unsigned i = 0; i < jets_handle->size(); ++i) {
    pat::Jet const& src = jets_handle->at(i);
    ac::PFJet& dest = output()->at(i);

    dest.setVector(setVar("p4", src.polarP4()));
    dest.setCharge(setVar("charge", src.charge()));

    dest.setPassesJetID(jetidSelector_(src));
  }
}

DEFINE_FWK_MODULE(AcornPFJetProducer);
