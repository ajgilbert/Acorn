#include "Acorn/NTupler/plugins/AcornPhotonProducer.h"
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
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

AcornPhotonProducer::AcornPhotonProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::Photon>>(config),
      inputToken_(consumes<edm::View<reco::Photon>>(config.getParameter<edm::InputTag>("input"))),
      phoLooseIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("phoLooseIdMap"))),
      phoMediumIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("phoMediumIdMap"))),
      phoTightIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("phoTightIdMap"))) {}
// phoMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult>>(
//     config.getParameter<edm::InputTag>("phoMediumIdFullInfoMap")))

AcornPhotonProducer::~AcornPhotonProducer() { ;}

void AcornPhotonProducer::produce(edm::Event& event,
                               const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::Photon> > photons_handle;

  event.getByToken(inputToken_, photons_handle);

  edm::Handle<edm::ValueMap<bool>> looseIdDecisions;
  edm::Handle<edm::ValueMap<bool>> mediumIdDecisions;
  edm::Handle<edm::ValueMap<bool>> tightIdDecisions;
  event.getByToken(phoLooseIdMapToken_, looseIdDecisions);
  event.getByToken(phoMediumIdMapToken_, mediumIdDecisions);
  event.getByToken(phoTightIdMapToken_, tightIdDecisions);
  // edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_data;
  // event.getByToken(phoMediumIdFullInfoMapToken_,medium_id_cutflow_data);

  output()->clear();
  output()->resize(photons_handle->size(), ac::Photon());

  for (unsigned i = 0; i < photons_handle->size(); ++i) {
    reco::Photon const& src = photons_handle->at(i);
    ac::Photon& dest = output()->at(i);

    dest.setIsLooseIdPhoton(
        setVar("isLooseIdPhoton", (*looseIdDecisions)[photons_handle->ptrAt(i)]));
    dest.setIsMediumIdPhoton(
        setVar("isMediumIdPhoton", (*mediumIdDecisions)[photons_handle->ptrAt(i)]));
    dest.setIsTightIdPhoton(
        setVar("isTightIdPhoton", (*tightIdDecisions)[photons_handle->ptrAt(i)]));

    // Two variables we need are only implemented for pat::Photons
    pat::Photon const* pat_src = dynamic_cast<pat::Photon const*>(&src);
    if (pat_src) {
      dest.setPassElectronVeto(setVar("passElectronVeto", pat_src->passElectronVeto()));
      dest.setHasPixelSeed(setVar("hasPixelSeed", pat_src->hasPixelSeed()));
    }


    // This is how to print the full cut flow later if desired.
    // See this message for details of what the cryptic cut names are:
    // https://hypernews.cern.ch/HyperNews/CMS/get/met/506/2.html
    // vid::CutFlowResult fullCutFlowData = (*medium_id_cutflow_data)[photons_handle->ptrAt(i)];
    // printCutFlowResult(fullCutFlowData);


    dest.setVector(setVar("p4", src.polarP4()));
    dest.setCharge(setVar("charge", src.charge()));
  }
}

DEFINE_FWK_MODULE(AcornPhotonProducer);
