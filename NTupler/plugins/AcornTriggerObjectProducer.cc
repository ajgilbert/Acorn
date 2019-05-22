#include "Acorn/NTupler/plugins/AcornTriggerObjectProducer.h"
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include "boost/format.hpp"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/NTupler/interface/city.h"

AcornTriggerObjectProducer::AcornTriggerObjectProducer(edm::ParameterSet const& config)
    : AcornBaseProducer<std::vector<ac::TriggerObject>>(config),
      inputToken_(consumes<pat::TriggerObjectStandAloneCollection>(
          config.getParameter<edm::InputTag>("input"))),
      triggerResultsToken_(
          consumes<edm::TriggerResults>(config.getParameter<edm::InputTag>("triggerResults"))),
      hltConfigProcess_(config.getParameter<std::string>("hltConfigProcess")),
      hltPath_(config.getParameter<std::string>("hltPath")),
      storeIfFired_(config.getParameter<bool>("storeIfFired")),
      extraFilterLabels_(config.getParameter<std::vector<std::string>>("extraFilterLabels")) {}

AcornTriggerObjectProducer::~AcornTriggerObjectProducer() { ; }

void AcornTriggerObjectProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  // We need to figure out the full HLT path name (typically hlt_path_ doesn't
  // contain the version number at the end). We don't have access to the full
  // pat::TriggerEvent here, but we can get the same information from the
  // edm::TriggerResults instead.
  // This code was inspired by the example here:
  //  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Trigger

  edm::Handle<pat::TriggerObjectStandAloneCollection> trigobj_handle;
  event.getByToken(inputToken_, trigobj_handle);

  edm::Handle<edm::TriggerResults> trigres_handle;
  event.getByToken(triggerResultsToken_, trigres_handle);
  edm::TriggerNames const& names = event.triggerNames(*trigres_handle);

  output()->clear();

  bool fired = true;
  bool path_found = false;
  std::string full_name;

  for (unsigned int i = 0, n = trigres_handle->size(); i < n; ++i) {
    std::string const& name = names.triggerName(i);
    // std::cout << i << "\t" << name << "\n";
    if (name.find(hltPath_) != name.npos) {
      full_name = name;
      path_found = true;
      if (storeIfFired_ && !(trigres_handle->accept(i))) fired = false;
      break;  // Stop loop after we find the first match
    }
  }
  if (!fired || !path_found) return;

  // Have to use the HLTConfigProvider to get the list of object-producing
  // filter modules that were run in this path
  std::set<std::string> path_filters;
  std::vector<std::string> const& filt_vec = hltConfig_.saveTagsModules(full_name);
  for (unsigned i = 0; i < filt_vec.size(); ++i) path_filters.insert(filt_vec[i]);
  for (unsigned i = 0; i < extraFilterLabels_.size(); ++i) path_filters.insert(extraFilterLabels_[i]);

  for (unsigned i = 0; i < trigobj_handle->size(); ++i) {
    pat::TriggerObjectStandAlone src = trigobj_handle->at(i);

    src.unpackPathNames(names);

    std::vector<std::string> const& pathnames = src.pathNames(false, false);
    bool obj_in_path = false;
    for (unsigned j = 0; j < pathnames.size(); ++j) {
      if (full_name == pathnames[j]) {
        obj_in_path = true;
        break;
      }
    }
    if (!obj_in_path) continue;

    // Ok, so this object was used in the path - now we can add it to the
    // output collection
    output()->push_back(ac::TriggerObject());
    ac::TriggerObject& dest = output()->back();
    dest.setVector(setVar("p4", src.polarP4()));

    // Get the filters this object was used in
// From 9_X_Y onwards have to unpack the labels too
#if CMSSW_MAJOR_VERSION >= 9
    src.unpackFilterLabels(event, *trigres_handle);
#endif
    std::vector<std::string> const& filters = src.filterLabels();
    std::vector<uint64_t> filter_labels;
    for (unsigned k = 0; k < filters.size(); ++k) {
      // Using the info we got from the HLTConfigProvider we can check if this
      // filter module was actually used in the path we are interested in
      if (!path_filters.count(filters[k])) continue;
      filter_labels.push_back(CityHash64(filters[k]));
      observedFilters_[filters[k]] = CityHash64(filters[k]);
    }
    dest.setFilters(filter_labels);

    // Assuming we can represent each trigger type as a short int, see here:
    // github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/HLTReco/interface/TriggerTypeDefs.h
    // we can pack four of these, each 16 bits, into the id() variable
    // (64 bits).
    ui64 packed_type;
    packed_type.one = 0;
    unsigned n_types = std::min(std::size_t(4), src.triggerObjectTypes().size());
    for (unsigned t = 0; t < n_types; ++t) {
      packed_type.four[t] = src.triggerObjectTypes()[t];
    }
    dest.setId(packed_type.one);
  }
}

void AcornTriggerObjectProducer::beginRun(edm::Run const& run, edm::EventSetup const& es) {
  AcornBaseProducer::beginRun(run, es);
  bool changed = true;
  bool res = hltConfig_.init(run, es, hltConfigProcess_, changed);
  if (!res) throw cms::Exception("HLTConfigProvider did not initialise correctly");
  observedHltMenus_.insert(hltConfig_.tableName());
}

void AcornTriggerObjectProducer::endStream() {
  std::cout << std::string(78, '-') << "\n";
  std::cout << boost::format("Path: %-50s  %20s\n")
      % hltPath_ % std::string("Hash Summmary");
  std::map<std::string, std::size_t>::const_iterator iter;
  for (iter = observedFilters_.begin(); iter != observedFilters_.end();
       ++iter) {
    // ICHashTreeProducer::Add(iter->second, iter->first);
    std::cout << boost::format("%-56s| %020i\n") % iter->first % iter->second;
  }
  std::cout << boost::format("HLT menus:\n");
  for (std::string const& str : observedHltMenus_) {
    std::cout << boost::format("%s\n") % str;
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornTriggerObjectProducer);
