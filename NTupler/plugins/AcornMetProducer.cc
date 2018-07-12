#include "Acorn/NTupler/plugins/AcornMetProducer.h"
#include <memory>
#include <string>
#include <vector>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/MET.h"
#include "Acorn/NTupler/interface/Met.h"

AcornMetProducer::AcornMetProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<std::vector<ac::Met>>(config),
      inputToken_(consumes<edm::View<reco::MET>>(config.getParameter<edm::InputTag>("input"))),
      saveGenMetFromPat_(config.getParameter<bool>("saveGenMetFromPat")) {}

AcornMetProducer::~AcornMetProducer() { ; }

void AcornMetProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::MET>> met_handle;
  event.getByToken(inputToken_, met_handle);

  output()->clear();
  output()->resize(met_handle->size(), ac::Met());

  for (unsigned i = 0; i < met_handle->size(); ++i) {
    reco::MET const& src = met_handle->at(i);

    // Below we will use a pointer to access the MET properties. This pointer might be changed to a
    // different reco::MET object first however, e.g. if the user wants to save the genMET that is
    // embedded in the pat::MET object
    reco::MET const* srcptr = &src;


    if (saveGenMetFromPat_) {
      pat::MET const* patsrc = dynamic_cast<pat::MET const*>(&src);
      if (!patsrc) {
        throw cms::Exception("DynamicCastFailed") << "Failed to dynamic cast reco::MET object to pat::MET\n";
      }
      srcptr = patsrc->genMET();
      if (!srcptr) {
        throw cms::Exception("NoGenMETAvailable") << "pat::MET object did not contain the genMET\n";
      }
    }

    ac::Met & dest = output()->at(i);

    dest.setVector(setVar("p4", srcptr->polarP4()));
    dest.setSumEt(setVar("sumEt", srcptr->sumEt()));
  }
}

DEFINE_FWK_MODULE(AcornMetProducer);
