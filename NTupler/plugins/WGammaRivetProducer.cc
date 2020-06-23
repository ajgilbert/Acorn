#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "Rivet/AnalysisHandler.hh"
#include "Acorn/NTupler/src/CMS_2020_PAS_SMP_20_005.h"

#include <vector>
#include <cstdio>
#include <cstring>

using namespace Rivet;
using namespace edm;

class WGammaRivetProducer : public edm::one::EDProducer<edm::one::SharedResources> {
public:
  explicit WGammaRivetProducer(const edm::ParameterSet& cfg)
      : _hepmcCollection(consumes<HepMCProduct>(cfg.getParameter<edm::InputTag>("HepMCCollection"))) {
    usesResource("Rivet");
    _analysis = new Rivet::CMS_2020_PAS_SMP_20_005();
    _isFirstEvent = true;
    produces<WGammaRivetVariables>();
  }

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<edm::HepMCProduct> _hepmcCollection;

  Rivet::AnalysisHandler _analysisHandler;
  Rivet::CMS_2020_PAS_SMP_20_005* _analysis;

  bool _isFirstEvent;
};

void WGammaRivetProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<HepMCProduct> evt;

  bool product_exists = iEvent.getByToken(_hepmcCollection, evt);
  if (product_exists) {
    HepMC::GenEvent const* myGenEvent = evt->GetEvent();

    if(_isFirstEvent) {
      _analysisHandler.addAnalysis(_analysis);
      _analysisHandler.init(*myGenEvent);
      _isFirstEvent = false;
    }
    _analysisHandler.analyze(*myGenEvent);
  }

  std::unique_ptr<WGammaRivetVariables> vars(new WGammaRivetVariables(_analysis->vars_));

  iEvent.put(std::move(vars));
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(WGammaRivetProducer);
