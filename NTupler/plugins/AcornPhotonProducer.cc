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
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("phoTightIdMap"))),
      phoCutFlowToken_(consumes<edm::ValueMap<vid::CutFlowResult>>(
          config.getParameter<edm::InputTag>("phoCutFlow"))),
      chargedIsolationLabel_(config.getParameter<std::string>("chargedIsolation")),
      neutralHadronIsolationLabel_(config.getParameter<std::string>("neutralHadronIsolation")),
      photonIsolationLabel_(config.getParameter<std::string>("photonIsolation")) {}

AcornPhotonProducer::~AcornPhotonProducer() { ;}

void AcornPhotonProducer::produce(edm::Event& event,
                               const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::Photon> > photons_handle;

  event.getByToken(inputToken_, photons_handle);

  edm::Handle<edm::ValueMap<bool>> looseIdDecisions;
  edm::Handle<edm::ValueMap<bool>> mediumIdDecisions;
  edm::Handle<edm::ValueMap<bool>> tightIdDecisions;
  edm::Handle<edm::ValueMap<vid::CutFlowResult>> cutflow;

  event.getByToken(phoLooseIdMapToken_, looseIdDecisions);
  event.getByToken(phoMediumIdMapToken_, mediumIdDecisions);
  event.getByToken(phoTightIdMapToken_, tightIdDecisions);
  event.getByToken(phoCutFlowToken_, cutflow);

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

    dest.setScEta(setVar("scEta", src.superCluster()->eta()));
    dest.setHadTowOverEm(setVar("hadTowOverEm", src.hadTowOverEm()));
    dest.setFull5x5SigmaIetaIeta(setVar("full5x5SigmaIetaIeta", src.full5x5_sigmaIetaIeta()));

    // This is how to print the full cut flow later if desired.
    // See this message for details of what the cryptic cut names are:
    // https://hypernews.cern.ch/HyperNews/CMS/get/met/506/2.html
    vid::CutFlowResult fullCutFlowData = (*cutflow)[photons_handle->ptrAt(i)];
    // printCutFlowResult(fullCutFlowData);

    dest.setChargedIso(setVar("chargedIso", fullCutFlowData.getValueCutUpon(chargedIsolationLabel_)));
    dest.setNeutralHadronIso(setVar("neutralHadronIso", fullCutFlowData.getValueCutUpon(neutralHadronIsolationLabel_)));
    dest.setPhotonIso(setVar("photonIso", fullCutFlowData.getValueCutUpon(photonIsolationLabel_)));

    dest.setVector(setVar("p4", src.polarP4()));
    dest.setCharge(setVar("charge", src.charge()));
  }
}

void AcornPhotonProducer::printCutFlowResult(vid::CutFlowResult& cutflow) {
  printf("    CutFlow name= %s    decision is %d\n", cutflow.cutFlowName().c_str(),
         (int)cutflow.cutFlowPassed());
  int ncuts = cutflow.cutFlowSize();
  printf(
      " Index                               cut name              isMasked    value-cut-upon     "
      "pass?\n");
  for (int icut = 0; icut < ncuts; icut++) {
    printf("  %d       %50s    %d        %f          %d\n", icut,
           cutflow.getNameAtIndex(icut).c_str(), (int)cutflow.isCutMasked(icut),
           cutflow.getValueCutUpon(icut), (int)cutflow.getCutResultByIndex(icut));
  }
}

DEFINE_FWK_MODULE(AcornPhotonProducer);
