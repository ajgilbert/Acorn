#ifndef Acorn_NTupler_AcornTrackProducer_h
#define Acorn_NTupler_AcornTrackProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/Track.h"

class AcornTrackProducer : public AcornBaseProducer<std::vector<ac::Track>> {
 public:
  explicit AcornTrackProducer(edm::ParameterSet const &config);
  ~AcornTrackProducer();

 private:
  virtual void produce(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<edm::View<reco::Track>> inputToken_;
  bool storeSlimmed_;
  boost::hash<reco::Track const*> track_hasher_;
};

#endif

