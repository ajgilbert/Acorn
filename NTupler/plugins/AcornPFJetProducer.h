#ifndef Acorn_NTupler_AcornPFJetProducer_h
#define Acorn_NTupler_AcornPFJetProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/PFJet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

class AcornPFJetProducer : public AcornBaseProducer<std::vector<ac::PFJet>> {
 public:
  explicit AcornPFJetProducer(edm::ParameterSet const &config);
  ~AcornPFJetProducer();

 private:
  virtual void produce(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<edm::View<pat::Jet>> inputToken_;

  PFJetIDSelectionFunctor jetidSelector_;
};

#endif
