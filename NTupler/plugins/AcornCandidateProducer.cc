#include "Acorn/NTupler/plugins/AcornCandidateProducer.h"
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
#include "Acorn/NTupler/interface/Candidate.h"

AcornCandidateProducer::AcornCandidateProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::Candidate>>(config),
      inputToken_(
          consumes<edm::View<reco::Candidate>>(config.getParameter<edm::InputTag>("input"))) {}

AcornCandidateProducer::~AcornCandidateProducer() { ; }

void AcornCandidateProducer::produce(edm::Event& event,
                                    const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::Candidate> > parts_handle;
  event.getByToken(inputToken_, parts_handle);

  output()->clear();
  output()->resize(parts_handle->size(), ac::Candidate());

  for (unsigned i = 0; i < parts_handle->size(); ++i) {
    reco::Candidate const& src = parts_handle->at(i);
    ac::Candidate& dest = output()->at(i);

    dest.setVector(setVar("p4", src.polarP4()));
    dest.setCharge(setVar("charge", src.charge()));
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornCandidateProducer);
