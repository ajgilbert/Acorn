#include "Acorn/NTupler/plugins/AcornMetProducer.h"
#include <memory>
#include <string>
#include <vector>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/MET.h"
#include "Acorn/NTupler/interface/Met.h"

AcornMetProducer::AcornMetProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::Met>>(config),
      inputToken_(
          consumes<edm::View<reco::MET>>(config.getParameter<edm::InputTag>("input"))) {}

AcornMetProducer::~AcornMetProducer() { ; }

void AcornMetProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::MET>> met_handle;
  event.getByToken(inputToken_, met_handle);

  output()->clear();
  output()->resize(met_handle->size(), ac::Met());

  for (unsigned i = 0; i < met_handle->size(); ++i) {
    reco::MET const& src = met_handle->at(i);
    ac::Met & dest = output()->at(i);

    dest.setVector(setVar("p4", src.polarP4()));
    dest.setSumEt(setVar("sumEt", src.sumEt()));
  }
}

DEFINE_FWK_MODULE(AcornMetProducer);
