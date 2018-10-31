#include "Acorn/NTupler/plugins/AcornElectronProducer.h"
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

AcornElectronProducer::AcornElectronProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::Electron>>(config),
      inputToken_(
          consumes<edm::View<reco::GsfElectron>>(config.getParameter<edm::InputTag>("input"))),
      vertexToken_(
          consumes<edm::View<reco::Vertex>>(config.getParameter<edm::InputTag>("inputVertices"))),
      eleVetoIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("eleVetoIdMap"))),
      eleLooseIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("eleLooseIdMap"))),
      eleMediumIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("eleMediumIdMap"))),
      eleTightIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("eleTightIdMap"))),
      eleMVAwp80IdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("eleMVAwp80IdMap"))),
      eleMVAwp90IdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("eleMVAwp90IdMap"))),
      eleHEEPIdMapToken_(
          consumes<edm::ValueMap<bool>>(config.getParameter<edm::InputTag>("eleHEEPIdMap"))) {}


AcornElectronProducer::~AcornElectronProducer() { ; }

void AcornElectronProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::GsfElectron>> electrons_handle;
  event.getByToken(inputToken_, electrons_handle);

  edm::Handle<edm::ValueMap<bool>> vetoIdDecisions;
  edm::Handle<edm::ValueMap<bool>> looseIdDecisions;
  edm::Handle<edm::ValueMap<bool>> mediumIdDecisions;
  edm::Handle<edm::ValueMap<bool>> tightIdDecisions;
  edm::Handle<edm::ValueMap<bool>> mvawp80IdDecisions;
  edm::Handle<edm::ValueMap<bool>> mvawp90IdDecisions;
  edm::Handle<edm::ValueMap<bool>> heepIdDecisions;

  event.getByToken(eleVetoIdMapToken_, vetoIdDecisions);
  event.getByToken(eleLooseIdMapToken_, looseIdDecisions);
  event.getByToken(eleMediumIdMapToken_, mediumIdDecisions);
  event.getByToken(eleTightIdMapToken_, tightIdDecisions);
  event.getByToken(eleMVAwp80IdMapToken_, mvawp80IdDecisions);
  event.getByToken(eleMVAwp90IdMapToken_, mvawp90IdDecisions);
  event.getByToken(eleHEEPIdMapToken_, heepIdDecisions);

  output()->clear();
  output()->resize(electrons_handle->size(), ac::Electron());

  edm::Handle<edm::View<reco::Vertex>> vtx_handle;
  event.getByToken(vertexToken_, vtx_handle);

  // It's possible we have no vertex, in which case give
  // it a dummy default one at (0, 0, 0)
  reco::Vertex dummyVertex;
  reco::Vertex const* firstVertex = &dummyVertex;
  if (vtx_handle->size() > 0) {
    firstVertex = &(vtx_handle->at(0));
  }

  for (unsigned i = 0; i < electrons_handle->size(); ++i) {
    reco::GsfElectron const& src = electrons_handle->at(i);
    ac::Electron& dest = output()->at(i);

    dest.setVector(setVar("p4", src.polarP4()));
    dest.setCharge(setVar("charge", src.charge()));

    dest.setIsCutBasedVetoElectron(setVar("isCutBasedVetoElectron",(*vetoIdDecisions)[electrons_handle->ptrAt(i)]));
    dest.setIsCutBasedLooseElectron(setVar("isCutBasedLooseElectron",(*looseIdDecisions)[electrons_handle->ptrAt(i)]));
    dest.setIsCutBasedMediumElectron(setVar("isCutBasedLooseElectron",(*mediumIdDecisions)[electrons_handle->ptrAt(i)]));
    dest.setIsCutBasedTightElectron(setVar("isCutBasedLooseElectron",(*tightIdDecisions)[electrons_handle->ptrAt(i)]));

    dest.setIsMVAwp80Electron(setVar("isMVAwp80Electron",(*mvawp80IdDecisions)[electrons_handle->ptrAt(i)]));
    dest.setIsMVAwp90Electron(setVar("isMVAwp90Electron",(*mvawp90IdDecisions)[electrons_handle->ptrAt(i)]));

    dest.setIsHEEPElectron(setVar("isHEEPElectron",(*heepIdDecisions)[electrons_handle->ptrAt(i)]));

    dest.setDxy(setVar("dxy", src.gsfTrack()->dxy(firstVertex->position())));
    dest.setDz(setVar("dz", src.gsfTrack()->dz(firstVertex->position())));

    dest.setVertex(setVar("vertex", src.vertex()));

  }
}

DEFINE_FWK_MODULE(AcornElectronProducer);
