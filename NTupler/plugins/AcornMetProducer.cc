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
      saveCorrectionLevels_(config.getParameter<std::vector<int>>("saveCorrectionLevels")),
      saveUncertaintyShifts_(config.getParameter<std::vector<int>>("saveUncertaintyShifts")),
      saveGenMetFromPat_(config.getParameter<bool>("saveGenMetFromPat")),
      skipMainMet_(config.getParameter<bool>("skipMainMet")) {}

AcornMetProducer::~AcornMetProducer() { ; }

void AcornMetProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::MET>> met_handle;
  event.getByToken(inputToken_, met_handle);

  // Number of MET objects to save from the main collection
  unsigned n_from_main = skipMainMet_ ? 0 : met_handle->size();

  output()->clear();
  output()->resize(n_from_main, ac::Met());

  for (unsigned i = 0; i < n_from_main; ++i) {
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

  if (saveCorrectionLevels_.size()) {
    pat::MET const* patsrc = dynamic_cast<pat::MET const*>(&met_handle->at(0));
    if (!patsrc) {
      throw cms::Exception("DynamicCastFailed")
          << "Failed to dynamic cast reco::MET object to pat::MET\n";
    }
    for (unsigned icorr = 0; icorr < saveCorrectionLevels_.size(); ++icorr) {
      auto level = pat::MET::METCorrectionLevel(saveCorrectionLevels_[icorr]);
      output()->push_back(ac::Met());
      ac::Met& dest = output()->back();
      dest.setLevel(level);
      dest.setVector(setVar("p4", ROOT::Math::PtEtaPhiMVector(patsrc->corPt(level), 0., patsrc->corPhi(level), 0.)));
      dest.setSumEt(setVar("sumEt", patsrc->corSumEt(level)));
      // std::cout << level << "\t" << dest.vector() << "\n";
      for (unsigned ishift = 0; ishift < saveUncertaintyShifts_.size(); ++ishift) {
        auto shift = pat::MET::METUncertainty(saveUncertaintyShifts_[ishift]);
        output()->push_back(ac::Met());
        ac::Met& dest = output()->back();
        dest.setLevel(level);
        dest.setShift(shift);
        dest.setVector(setVar("p4", ROOT::Math::PtEtaPhiMVector(patsrc->shiftedPt(shift, level), 0., patsrc->shiftedPhi(shift, level), 0.)));
        dest.setSumEt(setVar("sumEt", patsrc->shiftedSumEt(shift, level)));
        // std::cout << level << "\t" << shift << "\t" << dest.vector() << "\n";
      }
    }
  }
}

DEFINE_FWK_MODULE(AcornMetProducer);
