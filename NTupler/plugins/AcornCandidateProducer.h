#ifndef Acorn_NTupler_AcornCandidateProducer_h
#define Acorn_NTupler_AcornCandidateProducer_h

#include <memory>
#include <vector>
#include <string>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Rtypes.h"
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/GenParticle.h"


class AcornCandidateProducer : public AcornBaseProducer<std::vector<ac::Candidate>> {
 public:
  explicit AcornCandidateProducer(edm::ParameterSet const &config);
  ~AcornCandidateProducer();

 private:
  virtual void produce(edm::Event &, const edm::EventSetup &);

  edm::EDGetTokenT<edm::View<reco::Candidate>> inputToken_;
};



#endif
