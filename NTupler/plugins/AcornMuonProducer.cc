#include "Acorn/NTupler/plugins/AcornMuonProducer.h"
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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "Acorn/NTupler/interface/Muon.h"

AcornMuonProducer::AcornMuonProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::Muon>>(config),
      inputToken_(
          consumes<edm::View<reco::Muon>>(config.getParameter<edm::InputTag>("input"))),
      vertexToken_(
          consumes<edm::View<reco::Vertex>>(config.getParameter<edm::InputTag>("inputVertices"))) {}

AcornMuonProducer::~AcornMuonProducer() { ; }

void AcornMuonProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::Muon>> muons_handle;
  event.getByToken(inputToken_, muons_handle);

  output()->clear();
  output()->resize(muons_handle->size(), ac::Muon());

  edm::Handle<edm::View<reco::Vertex>> vtx_handle;
  event.getByToken(vertexToken_, vtx_handle);

  // It's possible we have no vertex, in which case give
  // it a dummy default one at (0, 0, 0)
  reco::Vertex dummyVertex;
  reco::Vertex const* firstVertex = &dummyVertex;
  if (vtx_handle->size() > 0) {
    firstVertex = &(vtx_handle->at(0));
  }

  for (unsigned i = 0; i < muons_handle->size(); ++i) {
    reco::Muon const& src = muons_handle->at(i);
    ac::Muon& dest = output()->at(i);

    dest.setVector(setVar("p4", src.polarP4()));
    dest.setCharge(setVar("charge", src.charge()));

    dest.setIsMediumMuon(setVar("isMediumMuon", muon::isMediumMuon(src)));
    dest.setIsTightMuon(setVar("isTightMuon", muon::isTightMuon(src, *firstVertex)));

    dest.setDxy(setVar("dxy", src.muonBestTrack()->dxy(firstVertex->position())));
    dest.setDz(setVar("dz", src.muonBestTrack()->dz(firstVertex->position())));

    dest.setVertex(setVar("vertex", src.vertex()));

    dest.SetPfIsoSumChargedHadronPt(
        setVar("pfIsoSumChargedHadronPt", src.pfIsolationR04().sumChargedHadronPt));
    dest.SetPfIsoSumNeutralHadronEt(
        setVar("pfIsoSumNeutralHadronEt", src.pfIsolationR04().sumNeutralHadronEt));
    dest.SetPfIsoSumPhotonEt(
        setVar("pfIsoSumPhotonEt", src.pfIsolationR04().sumPhotonEt));
    dest.SetPfIsoSumPUPt(
        setVar("pfIsoSumPUPt", src.pfIsolationR04().sumPUPt));
  }
}

DEFINE_FWK_MODULE(AcornMuonProducer);
