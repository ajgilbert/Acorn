#include "Acorn/NTupler/plugins/ICGenParticleProducer.hh"
#include <string>
#include <vector>
#include <bitset>
#include "TLorentzVector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

ICGenParticleProducer::ICGenParticleProducer(const edm::ParameterSet& config)
    : input_(config.getParameter<edm::InputTag>("input")),
      branch_(config.getParameter<std::string>("branch")),
      store_mothers_(config.getParameter<bool>("includeMothers")),
      store_daughters_(config.getParameter<bool>("includeDaughters")),
      store_statusFlags_(config.getParameter<bool>("includeStatusFlags")){
  consumes<edm::View<reco::GenParticle>>(input_);
  particles_ = new std::vector<ROOT::Math::PtEtaPhiMVector>();

  varActions_["p4"] = {false, true, -1};

  // PrintHeaderWithProduces(config, input_, branch_);
  // PrintOptional(1, store_mothers_, "includeMothers");
  // PrintOptional(1, store_daughters_, "includeDaughters");
  // PrintOptional(1, store_statusFlags_, "includeStatusFlags");
}

ICGenParticleProducer::~ICGenParticleProducer() { delete particles_; }

void ICGenParticleProducer::produce(edm::Event& event,
                                    const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::GenParticle> > parts_handle;
  event.getByLabel(input_, parts_handle);

  particles_->clear();
  particles_->resize(parts_handle->size(), ROOT::Math::PtEtaPhiMVector());

  for (unsigned i = 0; i < parts_handle->size(); ++i) {
    reco::GenParticle const& src = parts_handle->at(i);
    // edm::RefToBase<reco::GenParticle> ref = parts_handle->refAt(i);
    ROOT::Math::PtEtaPhiMVector& dest = particles_->at(i);
    // dest.SetPt(reduceMantissaToNbitsD(src.pt(), 12));
    // dest.SetEta(reduceMantissaToNbitsD(src.eta(), 12));
    // dest.SetPhi(reduceMantissaToNbitsD(src.phi(), 12));
    // dest.SetM(reduceMantissaToNbitsD(src.mass(), 12));

    dest = SetVar("p4", src.polarP4());
    // dest.set_id(particle_hasher_(&src));
    // dest.setPt(src.pt());
    // dest.setEta(src.eta());
    // dest.setPhi(src.phi());
    // dest.setE(src.energy());
    // dest.SetPtEtaPhiE(src.pt(), src.eta(), src.phi(), src.energy());
  }
}


// define this as a plug-in
DEFINE_FWK_MODULE(ICGenParticleProducer);
