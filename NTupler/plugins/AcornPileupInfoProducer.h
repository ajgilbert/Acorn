#ifndef Acorn_NTupler_AcornPileupInfoProducer_h
#define Acorn_NTupler_AcornPileupInfoProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/PileupInfo.h"

class AcornPileupInfoProducer : public AcornBaseProducer<std::vector<ac::PileupInfo>> {
 public:
  explicit AcornPileupInfoProducer(edm::ParameterSet const &config);
  ~AcornPileupInfoProducer();

 private:
  virtual void produce(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<edm::View<PileupSummaryInfo>> inputToken_;
  int minBx_;
  int maxBx_;
};

#endif
