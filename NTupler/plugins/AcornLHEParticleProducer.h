#ifndef Acorn_NTupler_AcornLHEParticleProducer_h
#define Acorn_NTupler_AcornLHEParticleProducer_h

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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

class AcornLHEParticleProducer : public AcornBaseProducer<std::vector<ac::GenParticle>> {
 public:
  explicit AcornLHEParticleProducer(const edm::ParameterSet &config);
  ~AcornLHEParticleProducer();

 private:
  virtual void produce(edm::Event &, const edm::EventSetup &);

  edm::EDGetTokenT<LHEEventProduct> inputToken_;
};



#endif
