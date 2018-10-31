#ifndef Acorn_NTupler_AcornElectronProducer_h
#define Acorn_NTupler_AcornElectronProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/Common/interface/ValueMap.h"

class AcornElectronProducer : public AcornBaseProducer<std::vector<ac::Electron>> {
 public:
  explicit AcornElectronProducer(edm::ParameterSet const &config);
  ~AcornElectronProducer();

 private:
  virtual void produce(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<edm::View<reco::GsfElectron>> inputToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> vertexToken_;

  edm::EDGetTokenT<edm::ValueMap<bool>> eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleTightIdMapToken_;

  edm::EDGetTokenT<edm::ValueMap<bool>> eleMVAwp80IdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleMVAwp90IdMapToken_;

  edm::EDGetTokenT<edm::ValueMap<bool>> eleHEEPIdMapToken_;

};

#endif
