#ifndef Acorn_NTupler_AcornMuonProducer_h
#define Acorn_NTupler_AcornMuonProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/Muon.h"

class AcornMuonProducer : public AcornBaseProducer<std::vector<ac::Muon>> {
 public:
  explicit AcornMuonProducer(edm::ParameterSet const &config);
  ~AcornMuonProducer();

 private:
  virtual void produce(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<edm::View<reco::Muon>> inputToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> vertexToken_;
};

#endif
