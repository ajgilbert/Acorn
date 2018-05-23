#include "Acorn/NTupler/plugins/AcornPileupInfoProducer.h"
#include <string>
#include <vector>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Acorn/NTupler/interface/PileupInfo.h"

AcornPileupInfoProducer::AcornPileupInfoProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::PileupInfo>>(config),
      inputToken_(
          consumes<edm::View<PileupSummaryInfo>>(config.getParameter<edm::InputTag>("input"))) {}

AcornPileupInfoProducer::~AcornPileupInfoProducer() { ; }

void AcornPileupInfoProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<PileupSummaryInfo>> info_handle;
  event.getByToken(inputToken_, info_handle);

  output()->clear();
  output()->resize(info_handle->size(), ac::PileupInfo());

  for (unsigned i = 0; i < info_handle->size(); ++i) {
    PileupSummaryInfo const& src = info_handle->at(i);
    ac::PileupInfo& dest = output()->at(i);

    dest.setNumInteractions(setVar("numInteractions", src.getPU_NumInteractions()));
    dest.setBunchCrossing(setVar("bunchCrossing", src.getBunchCrossing()));
    dest.setTrueNumInteractions(setVar("trueNumInteractions", src.getTrueNumInteractions()));
  }
}

DEFINE_FWK_MODULE(AcornPileupInfoProducer);
