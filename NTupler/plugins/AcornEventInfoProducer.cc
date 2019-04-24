#include "Acorn/NTupler/plugins/AcornEventInfoProducer.h"
#include <memory>
#include <string>
#include <vector>
#include <regex>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/StringUtils.h"
#include "FWCore/Utilities/interface/Exception.h"

AcornEventInfoProducer::AcornEventInfoProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<ac::EventInfo>(config),
      lheTag_(config.getParameter<edm::InputTag>("lheProducer")),
      lheToken_(consumes<LHEEventProduct>(config.getParameter<edm::InputTag>("lheProducer"))),
      genToken_(consumes<GenEventInfoProduct>(config.getParameter<edm::InputTag>("generator"))),
      includeLHEWeights_(config.getParameter<bool>("includeLHEWeights")),
      includeGenWeights_(config.getParameter<bool>("includeGenWeights")),
      metfilterToken_(consumes<edm::TriggerResults>(config.getParameter<edm::InputTag>("metFilterResults"))),
      saveMetFilters_(config.getParameter<std::vector<std::string>>("saveMetFilters")),
      includeNumVertices_(config.getParameter<bool>("includeNumVertices")),
      vertexToken_(consumes<edm::View<reco::Vertex>>(config.getParameter<edm::InputTag>("inputVertices"))) {
  consumes<LHERunInfoProduct, edm::InRun>({lheTag_});
  consumes<GenEventInfoProduct>({"generator"});

  for (auto const& tag : config.getParameter<std::vector<edm::InputTag>>("saveMetFilterBools")) {
    saveMetFilterBools_.emplace_back(consumes<bool>(tag));
  }

  for (auto const& tag : config.getParameter<std::vector<edm::InputTag>>("userDoubles")) {
    userDoubleTokens_.emplace_back(consumes<double>(tag));
  }
}

AcornEventInfoProducer::~AcornEventInfoProducer() {
  // delete info_;
}

void AcornEventInfoProducer::beginRun(edm::Run const & run, edm::EventSetup const &es) {
  AcornBaseProducer<ac::EventInfo>::beginRun(run, es);
  if (!includeLHEWeights_) return;
  if (lheWeightLabels_.size()) return;
  edm::Handle<LHERunInfoProduct> lhe_info;
  run.getByLabel(lheTag_, lhe_info);
  int record = 0;
  bool keepGroup = false;
  VarRule groupVarRule;
  // We need to be able to match something like:
  // <weight id="rwgt_100004">set param_card dim6 1 1.6 # orig: 0.0\n</weight>
  // Also need to catch id='X' (MG2.6.X) and id="X"
  std::regex rgx(R"(<weight.*id=[\"\'][^\d]*(\d+)[\"\'].*>([\s\S]*)</weight>)");
  std::regex rgx_group(R"(<weightgroup(.*)>)");
  for (auto it = lhe_info->headers_begin(); it != lhe_info->headers_end();
       ++it) {
    std::vector<std::string>::const_iterator iLt = it->begin();
    for (; iLt != it->end(); ++iLt) {
      std::string line = *iLt;
      // std::cout << line;
      // Fix for some headers produced with MG 2.6.X, the < and >
      // have been replaced with &lt; and &gt; everywhere
      boost::replace_all(line, "&lt;", "<");
      boost::replace_all(line, "&gt;", ">");
      if (line.find("<weightgroup")  != std::string::npos) {
        lheWeightLabels_.push_back(line);
        lheWeightWasKept_.push_back(false);
        // Have seen in MG2.6.X that a new weight group can start before the closing tag
        // therefore record is now an int, incremented for each <weightgroup> and
        // decremented for each </weightgroup>.
        record += 1;
        std::smatch rgx_match;
        std::regex_search(line, rgx_match, rgx_group);
        if (rgx_match.size() == 2) {
          // std::cout << " -- " << std::string(rgx_match[1]) << "\n";
          groupVarRule = getVarRule("lheweightgroups:" + std::string(rgx_match[1]));
          if (!groupVarRule.zeroed) {
            keepGroup = true;
            // std::cout << " -- Keeping group\n";
          }
        }
        continue;
      }
      if (line.find("</weightgroup") != std::string::npos) {
        lheWeightLabels_.push_back(line);
        lheWeightWasKept_.push_back(false);
        keepGroup = false;
        if (record > 0) record -= 1;
        continue;
      }
      if (record) {
        // For some entries madgraph adds a line break before the closing
        // </weight>, this is a workaround
        if (line.find("</weight>") == line.npos) {
          line = line + "</weight>";
        }
        // If this line contains "</weight>" at the beginning, it's just
        // overflow from the last line - we can skip it
        if (line.find("</weight>") == 0) {
          continue;
        }
        lheWeightLabels_.push_back(line);
        std::smatch rgx_match;
        std::regex_search(line, rgx_match, rgx);
        if (rgx_match.size() == 3) {
          std::vector<std::string> split_label = ac::TrimAndSplitString(rgx_match[2]);
          split_label.push_back(ac::TrimString(rgx_match[2])); // Also add the full string
          unsigned id = boost::lexical_cast<int>(rgx_match[1]);
          bool keep = false;
          if (keepGroup) {
            keep = true;
            savedLHEWeightIds[id] = groupVarRule;
          }
          for (auto const& x : split_label) {
            auto rule = getVarRule("lheweights:" + x);
            if (!rule.zeroed) {
              // std::cout << x << " matched " << rule.original << "\n";
              savedLHEWeightIds[id] = rule;
              keep = true;
              break;
            }
          }
          lheWeightWasKept_.push_back(keep);
        } else {
          lheWeightWasKept_.push_back(false);
          // Seen cases for madgraph 2.6.X where comment lines appear instead a <weightgroup>
          // block. If the string doesn't contain "<weight" we will just silently skip it
          if (line.find("<weight") == line.npos) {
            continue;
          }
          edm::LogWarning("LHEHeaderParsing") << "Line was not in the expected format: " << line << "\n";
        }
      }
    }
  }
}

void AcornEventInfoProducer::endRun(edm::Run const& run, edm::EventSetup const& es) {
}

void AcornEventInfoProducer::produce(edm::Event& event,
                                  const edm::EventSetup& setup) {
  ac::EventInfo * info = output();
  info->setIsRealData(event.isRealData());
  info->setRun(event.run());
  info->setEvent(event.id().event());
  info->setLuminosityBlock(event.luminosityBlock());
  info->setBunchCrossing(event.bunchCrossing());

  if (includeLHEWeights_) {
    edm::Handle<LHEEventProduct> lhe_handle;
    event.getByToken(lheToken_, lhe_handle);

    info->setNpLO(setVar("npLO", lhe_handle->npLO()));
    info->setNpNLO(setVar("npNLO", lhe_handle->npNLO()));

    double nominalLHEWeight = 1.;
    if (lhe_handle->weights().size()) {
      // Very rarely, this weights() vector is empty (seen in a powheg sample)
      nominalLHEWeight = lhe_handle->weights()[0].wgt;
    }
    info->setNominalLHEWeight(setVar("nominalLHEWeight", nominalLHEWeight));
    for (unsigned i = 0; i < lhe_handle->weights().size(); ++i) {
      // Weight id is a string, assume it can always cast to an unsigned int
      unsigned id = boost::lexical_cast<unsigned>(lhe_handle->weights()[i].id);
      auto it = savedLHEWeightIds.find(id);
      if (it != savedLHEWeightIds.end()) {
        double weight = lhe_handle->weights()[i].wgt / nominalLHEWeight;
        info->setLHEWeight(id, processVar(weight, it->second));
      }
    }
  }

  if (includeGenWeights_) {
    edm::Handle<GenEventInfoProduct> gen_info_handle;
    event.getByToken(genToken_, gen_info_handle);
    if (gen_info_handle.isValid()) {
      double nominal_gen_weight = gen_info_handle->weight();
      info->setNominalGenWeight(setVar("nominalGenWeight", gen_info_handle->weight()));

      // Insert the sign of the weight into the main weight calculator
      info->setWeight("wt_mc_sign", (nominal_gen_weight >= 0.) ? 1.0 : -1.0);

      std::vector<double> gen_weights(gen_info_handle->weights().size());
      for (unsigned i = 0; i < gen_info_handle->weights().size(); ++i) {
        gen_weights[i] = setVar("genWeights", gen_info_handle->weights()[i] / nominal_gen_weight);
      }
      info->setGenWeights(gen_weights);
    }
  }

  unsigned n_met_filters_to_save = saveMetFilters_.size() + saveMetFilterBools_.size();
  if (n_met_filters_to_save > 0) {
    constexpr unsigned maxfilters = 32;
    if (n_met_filters_to_save > maxfilters) {
      throw cms::Exception("TooManyMETFilters")
          << "Maximum MET filter flags to be saved is " << maxfilters << ", but "
          << n_met_filters_to_save << " were requested\n";
    }
    std::bitset<maxfilters> metfilter_bits;
    edm::Handle<edm::TriggerResults> metfilter_handle;
    event.getByToken(metfilterToken_, metfilter_handle);
    edm::TriggerNames const& triggerNames = event.triggerNames(*metfilter_handle);
    for (unsigned imet = 0; imet < saveMetFilters_.size(); ++imet) {
      auto trg_idx = triggerNames.triggerIndex(saveMetFilters_[imet]);
      if (trg_idx == triggerNames.size()) {
        // The string requested wasn't found, so we'll just skip it
        missedMetFilters_.insert(saveMetFilters_[imet]);
        continue;
      }
      // std::cout << imet << "\t" << triggerNames.triggerName(trg_idx) << "\t" << metfilter_handle->accept(trg_idx) << "\n";
      metfilter_bits[imet] = !metfilter_handle->accept(trg_idx);
    }
    for (unsigned imet = 0; imet < saveMetFilterBools_.size(); ++ imet) {
      edm::Handle<bool> bool_handle;
      event.getByToken(saveMetFilterBools_[imet], bool_handle);
      // Offset the index by the size of the metfilters saved from the TriggerResults
      metfilter_bits[saveMetFilters_.size() + imet] = !(*bool_handle);
    }
    info->setMetFilters(metfilter_bits);
    // std::cout << toBinaryString(uint32_t(metfilter_bits.to_ulong())) << "\n";
  }

  std::vector<double> user_doubles;
  for (unsigned idouble = 0; idouble < userDoubleTokens_.size(); ++idouble) {
    edm::Handle<double> double_handle;
    event.getByToken(userDoubleTokens_[idouble], double_handle);
    user_doubles.push_back(*double_handle);
  }
  info->setUserDoubles(user_doubles);

  if (includeNumVertices_) {
    edm::Handle<edm::View<reco::Vertex>> vtx_handle;
    event.getByToken(vertexToken_, vtx_handle);
    info->setNumVertices(vtx_handle->size());
  }
}

void AcornEventInfoProducer::endStream() {
  if (lheWeightLabels_.size()) {
    std::cout << "LHE event weights\n";
    std::cout << std::string(78, '-') << "\n";
    for (unsigned l = 0; l < lheWeightLabels_.size(); ++l) {
      if (lheWeightWasKept_[l]) std::cout << "(*) ";
      std::cout << lheWeightLabels_[l];
    }
  }
  if (missedMetFilters_.size()) {
    std::cout << "WARNING: the following MET filters were not defined for all events:\n";
    for (auto const& filtername : missedMetFilters_) {
      std::cout << " - " << filtername << "\n";
    }
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornEventInfoProducer);
