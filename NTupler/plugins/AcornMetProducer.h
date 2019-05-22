#ifndef Acorn_NTupler_AcornMetProducer_h
#define Acorn_NTupler_AcornMetProducer_h

#include <memory>
#include <vector>
#include <string>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"

class AcornMetProducer : public AcornBaseProducer<std::vector<ac::Met>> {
 public:
  explicit AcornMetProducer(edm::ParameterSet const &config);
  ~AcornMetProducer();

 private:
  virtual void produce(edm::Event&, edm::EventSetup const&);

  edm::EDGetTokenT<edm::View<reco::MET>> inputToken_;

  // Possible values:
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/DataFormats/PatCandidates/interface/MET.h#L151
  std::vector<int> saveCorrectionLevels_;
  std::vector<int> saveUncertaintyShifts_;
  bool saveGenMetFromPat_;
  bool skipMainMet_;
};


#endif
