#include "Acorn/NTupler/plugins/AcornRequestByDeltaR.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TrackReco/interface/Track.h"

typedef RequestByDeltaR<reco::Track,reco::Track> RequestTracksByDeltaRFromTrack;

DEFINE_FWK_MODULE(RequestTracksByDeltaRFromTrack);
