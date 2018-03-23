#ifndef Acorn_NTupler_AcornEventInfoProducer_h
#define Acorn_NTupler_AcornEventInfoProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/plugins/AcornBaseProducer.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

/**
 * @brief Produces an ac::EventInfo object
 *
 * **Example usage**
 * @snippet python/default_producers_cfi.py EventInfo
 */
class AcornEventInfoProducer : public AcornBaseProducer<ac::EventInfo> {
 public:
  explicit AcornEventInfoProducer(const edm::ParameterSet &);
  ~AcornEventInfoProducer();

 private:
  // virtual void beginJob();
  virtual void beginRun(edm::Run const& run, edm::EventSetup const& es);
  virtual void endRun(edm::Run const& run, edm::EventSetup const& es);
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endStream();


  std::map<unsigned, VarRule> savedLHEWeightIds;

  // ac::EventInfo *info_;
  // std::string branch_;
  edm::InputTag lheTag_;
  edm::EDGetTokenT<LHEEventProduct> lheToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genToken_;
  // bool do_jets_rho_;
  // edm::InputTag input_jets_rho_;
  // bool do_leptons_rho_;
  // edm::InputTag input_leptons_rho_;
  // bool do_vertex_count_;
  // edm::InputTag input_vertices_;
  bool includeLHEWeights_;
  bool includeGenWeights_;
  // bool do_embedding_weights_;
  // bool do_ht_;
  // std::vector<std::pair<std::string, edm::InputTag> > weights_;
  // std::vector<std::pair<std::string, edm::InputTag> > gen_weights_;
  std::vector<std::string> lheWeightLabels_;
  std::vector<bool> lheWeightWasKept_;

};

#endif
