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


  std::map<unsigned, VarRule> saved_lhe_weight_ids_;

  // ac::EventInfo *info_;
  // std::string branch_;
  edm::InputTag lhe_collection_;
  // bool do_jets_rho_;
  // edm::InputTag input_jets_rho_;
  // bool do_leptons_rho_;
  // edm::InputTag input_leptons_rho_;
  // bool do_vertex_count_;
  // edm::InputTag input_vertices_;
  bool do_lhe_weights_;
  // bool do_embedding_weights_;
  // bool do_ht_;
  // std::vector<std::pair<std::string, edm::InputTag> > weights_;
  // std::vector<std::pair<std::string, edm::InputTag> > gen_weights_;
  std::vector<std::string> lhe_weight_labels_;
  std::vector<bool> lhe_weight_was_kept_;

  //3 ways to get event filters:
  //1- CSC has special handle. Not active anymore for run II....
  // bool do_csc_filter_;
  // edm::InputTag input_csc_filter_;
  // //2- Bool filters from triggerResults handle, from miniAOD
  // bool do_filtersfromtrig_;
  // edm::InputTag filtersfromtrig_input_;
  // std::vector<std::string> filtersfromtrig_;
  // //3- Bool filters from specific collections - name starting with! will be inverted.
  // std::vector<std::pair<std::string, edm::InputTag> > filters_;
  // std::set<std::string> invert_filter_logic_;


  // //store all filters in one map
  // std::map<std::string, std::size_t> observed_filters_;
};

#endif
