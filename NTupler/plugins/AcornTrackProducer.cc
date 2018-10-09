#include "Acorn/NTupler/plugins/AcornTrackProducer.h"
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Acorn/NTupler/interface/Track.h"

AcornTrackProducer::AcornTrackProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::Track>>(config),
      inputToken_(consumes<edm::View<reco::Track>>(config.getParameter<edm::InputTag>("input"))){}

AcornTrackProducer::~AcornTrackProducer() { ; }

void AcornTrackProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::Track>> tracks_handle;
  event.getByToken(inputToken_, tracks_handle);

  output()->clear();
  output()->resize(tracks_handle->size(), ac::Track());

  for (unsigned i = 0; i < tracks_handle->size(); ++i) {
    reco::Track const& src = tracks_handle->at(i);
    ac::Track& dest = output()->at(i);

    dest.setPt(setVar("pt",src.pt()));
    dest.setEta(setVar("eta",src.eta()));
    dest.setPhi(setVar("phi",src.phi()));
    dest.setCharge(setVar("charge",src.charge()));
    dest.setVx(setVar("vx",src.vx()));
    dest.setVy(setVar("vy",src.vy()));
    dest.setVz(setVar("vz",src.vz()));
    dest.setHits(setVar("hits",src.hitPattern().numberOfValidHits()));
    dest.setPixel_hits(setVar("pixelhits",src.hitPattern().numberOfValidPixelHits()));
    dest.setQuality(setVar("qual",src.qualityMask()));
    #if CMSSW_MAJOR_VERSION >= 9
    dest.setHits_miss_inner(setVar("missinghits",
      src.hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)));
    #else 
    dest.setHits_miss_inner(setVar("missinghits",
      src.hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)));
    #endif
  }
}


DEFINE_FWK_MODULE(AcornTrackProducer);

