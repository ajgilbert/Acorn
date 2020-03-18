#ifndef Acorn_NTupler_AcornPhotonProducer_h
#define Acorn_NTupler_AcornPhotonProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/Common/interface/ValueMap.h"

class AcornPhotonProducer : public AcornBaseProducer<std::vector<ac::Photon>> {
 public:
  explicit AcornPhotonProducer(const edm::ParameterSet &);
  ~AcornPhotonProducer();

 private:
  virtual void produce(edm::Event &, const edm::EventSetup &);

  edm::EDGetTokenT<edm::View<reco::Photon>> inputToken_;

  edm::EDGetTokenT<edm::ValueMap<bool>> phoLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> phoMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> phoTightIdMapToken_;

  std::string phoLooseIdLabel_;
  std::string phoMediumIdLabel_;
  std::string phoTightIdLabel_;

  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > phoCutFlowToken_;

  std::string chargedIsolationLabel_;
  std::string neutralHadronIsolationLabel_;
  std::string photonIsolationLabel_;
  std::string worstChargedIsolationLabel_;

  // Take the id values stored within the pat::Photon objects
  // instead of reading from external ValueMaps
  bool takeIdsFromObjects_;

  std::vector<std::string> energyCorrectionLabels_;

  void printCutFlowResult(vid::CutFlowResult &cutflow);
};

#endif
