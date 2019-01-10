#ifndef Acorn_NTupler_AcornEventInfoProducer_h
#define Acorn_NTupler_AcornEventInfoProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

class AcornEventInfoProducer : public AcornBaseProducer<ac::EventInfo> {
 public:
  explicit AcornEventInfoProducer(const edm::ParameterSet &);
  ~AcornEventInfoProducer();

 private:
  // virtual void beginJob();
  virtual void beginRun(edm::Run const& run, edm::EventSetup const& es);
  virtual void endRun(edm::Run const& run, edm::EventSetup const& es);
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endStream();

  std::map<unsigned, VarRule> savedLHEWeightIds;

  edm::InputTag lheTag_;
  edm::EDGetTokenT<LHEEventProduct> lheToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genToken_;
  bool includeLHEWeights_;
  bool includeGenWeights_;
  std::vector<std::string> lheWeightLabels_;
  std::vector<bool> lheWeightWasKept_;
  edm::EDGetTokenT<edm::TriggerResults> metfilterToken_;
  std::vector<std::string> saveMetFilters_;
  std::set<std::string> missedMetFilters_;
  std::vector<edm::EDGetTokenT<double>> userDoubleTokens_;

};

#endif
