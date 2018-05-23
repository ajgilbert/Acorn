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

  // edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > phoMediumIdFullInfoMapToken_;
  // void printCutFlowResult(vid::CutFlowResult &cutflow) {
  //   printf("    CutFlow name= %s    decision is %d\n", cutflow.cutFlowName().c_str(),
  //          (int)cutflow.cutFlowPassed());
  //   int ncuts = cutflow.cutFlowSize();
  //   printf(
  //       " Index                               cut name              isMasked    value-cut-upon     "
  //       "pass?\n");
  //   for (int icut = 0; icut < ncuts; icut++) {
  //     printf("  %d       %50s    %d        %f          %d\n", icut,
  //            cutflow.getNameAtIndex(icut).c_str(), (int)cutflow.isCutMasked(icut),
  //            cutflow.getValueCutUpon(icut), (int)cutflow.getCutResultByIndex(icut));
  //   }
  // }
};

#endif
