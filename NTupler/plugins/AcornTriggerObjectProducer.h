#ifndef Acorn_NTupler_AcornObjectProducer_h
#define Acorn_NTupler_AcornObjectProducer_h

#include <memory>
#include <vector>
#include <set>
#include <string>
#include <cstdint>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"

#include "Acorn/NTupler/interface/TriggerObject.h"

class AcornTriggerObjectProducer : public AcornBaseProducer<std::vector<ac::TriggerObject>> {
 public:
  explicit AcornTriggerObjectProducer(const edm::ParameterSet&);
  ~AcornTriggerObjectProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginRun(edm::Run const& run, edm::EventSetup const& es);
  virtual void endStream();

  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> inputToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  std::string hltConfigProcess_;

  std::string hltPath_;
  bool storeIfFired_;
  std::map<std::string, std::size_t> observedFilters_;
  HLTConfigProvider hltConfig_;
  std::set<std::string> observedHltMenus_;
  std::vector<std::string> extraFilterLabels_;

  union ui64 {
    uint64_t one;
    int16_t four[4];
  };
};

#endif
