#include "Acorn/NTupler/plugins/AcornWGammaRivetProducer.h"
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

AcornWGammaRivetProducer::AcornWGammaRivetProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<WGammaRivetVariables>(config),
      inputToken_(
          consumes<WGammaRivetVariables>(config.getParameter<edm::InputTag>("input"))) {}

AcornWGammaRivetProducer::~AcornWGammaRivetProducer() { ; }

void AcornWGammaRivetProducer::produce(edm::Event& event,
                                    const edm::EventSetup& setup) {
  edm::Handle<WGammaRivetVariables> handle;
  event.getByToken(inputToken_, handle);

  *(output()) = *handle;
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornWGammaRivetProducer);
