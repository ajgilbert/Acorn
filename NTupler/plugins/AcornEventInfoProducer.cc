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
// #include "UserCode/ICHiggsTauTau/interface/StaticTree.hh"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/StringUtils.h"
// #include "Acorn/NTupler/interface/city.h"
// #include "Acorn/NTupler/plugins/PrintConfigTools.h"
#include "FWCore/Utilities/interface/Exception.h"
// #include "UserCode/ICHiggsTauTau/plugins/Consumes.h"

AcornEventInfoProducer::AcornEventInfoProducer(const edm::ParameterSet& config)
    : AcornBaseProducer<ac::EventInfo>(config),
      lheTag_(config.getParameter<edm::InputTag>("lheProducer")),
      lheToken_(consumes<LHEEventProduct>(config.getParameter<edm::InputTag>("lheProducer"))),
      genToken_(consumes<GenEventInfoProduct>(config.getParameter<edm::InputTag>("generator"))),
      includeLHEWeights_(config.getParameter<bool>("includeLHEWeights")),
      includeGenWeights_(config.getParameter<bool>("includeGenWeights"))
      // do_jets_rho_(config.getParameter<bool>("includeJetRho")),
      // input_jets_rho_(config.getParameter<edm::InputTag>("inputJetRho")),
      // do_leptons_rho_(config.getParameter<bool>("includeLeptonRho")),
      // input_leptons_rho_(config.getParameter<edm::InputTag>("inputLeptonRho")),
      // do_vertex_count_(config.getParameter<bool>("includeVertexCount")),
      // input_vertices_(config.getParameter<edm::InputTag>("inputVertices")),
// do_embedding_weights_(config.getParameter<bool>("includeEmbeddingWeights")),
// do_ht_(config.getParameter<bool>("includeHT")),
// do_csc_filter_(config.getParameter<bool>("includeCSCFilter")),
// input_csc_filter_(config.getParameter<edm::InputTag>("inputCSCFilter")),
// do_filtersfromtrig_(config.getParameter<bool>("includeFiltersFromTrig")),
// filtersfromtrig_input_(config.getParameter<edm::InputTag>("inputfiltersfromtrig")),
// filtersfromtrig_(config.getParameter<std::vector<std::string> >("filtersfromtrig"))
{
  consumes<LHERunInfoProduct, edm::InRun>({lheTag_});
  consumes<GenEventInfoProduct>({"generator"});
  // consumes<double>(input_leptons_rho_);
  // consumes<double>(input_jets_rho_);
  // consumes<edm::View<reco::Vertex>>(input_vertices_);
  // consumes<reco::BeamHaloSummary>(input_csc_filter_);
  // consumes<edm::TriggerResults>(filtersfromtrig_input_);
  // edm::ParameterSet filter_params =
  //     config.getParameter<edm::ParameterSet>("filters");
  // std::vector<std::string> filter_names =
  //     filter_params.getParameterNamesForType<edm::InputTag>();
  // for (unsigned i = 0; i < filter_names.size(); ++i) {
  //   filters_.push_back(std::make_pair(
  //       filter_names[i],
  //       filter_params.getParameter<edm::InputTag>(filter_names[i])));
  //       consumes<bool>(filters_[i].second);
  //   if (filters_.back().second.label().at(0) == '!') {
  //     std::cout << "Info in <AcornEventInfoProducer>: Inverting logic for filter: "
  //               << filters_.back().first << std::endl;
  //     std::string new_label = filters_.back().second.label();
  //     new_label.erase(0, 1);
  //     filters_.back().second =
  //         edm::InputTag(new_label, filters_.back().second.instance(),
  //                       filters_.back().second.process());
  //     invert_filter_logic_.insert(filters_.back().first);
  //   }
  // }
  // for (unsigned iFilter=0;iFilter<filtersfromtrig_.size();iFilter++){
  //   if((filtersfromtrig_[iFilter]).at(0)=='!'){
  //     filtersfromtrig_[iFilter].erase(0,1);
  //     invert_filter_logic_.insert(filtersfromtrig_[iFilter]);
  //   }
  // }

  // edm::ParameterSet wt_pset = config.getParameter<edm::ParameterSet>("weights");
  // std::vector<std::string> wt =
  //     wt_pset.getParameterNamesForType<edm::InputTag>();
  // for (unsigned i = 0; i < wt.size(); ++i) {
  //   weights_.push_back(
  //       std::make_pair(wt[i], wt_pset.getParameter<edm::InputTag>(wt[i])));
  //       consumes<double>(weights_[i].second);
  // }

  // edm::ParameterSet gwt_pset =
  //     config.getParameter<edm::ParameterSet>("genFilterWeights");
  // std::vector<std::string> gwt =
  //     gwt_pset.getParameterNamesForType<edm::InputTag>();
  // for (unsigned i = 0; i < gwt.size(); ++i) {
  //   gen_weights_.push_back(
  //       std::make_pair(gwt[i], gwt_pset.getParameter<edm::InputTag>(gwt[i])));
  //       consumes<double>(gen_weights_[i].second);
  // }

  // info_ = new ac::EventInfo();

  // PrintHeaderWithBranch(config, branch_);
  // PrintOptional(1, do_lhe_weights_, "includeLHEWeights");
  // PrintOptional(1, do_embedding_weights_, "includeEmbeddingWeights");
  // PrintOptional(1, do_ht_, "includeHT");
  // PrintOptional(1, do_jets_rho_, "includeJetRho");
  // PrintOptional(1, do_leptons_rho_, "includeLeptonRho");
  // PrintOptional(1, do_vertex_count_, "includeVertexCount");
  // PrintOptional(1, do_csc_filter_, "includeCSCFilter");
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
  // We need to be able to match something like:
  // <weight id="rwgt_100004">set param_card dim6 1 1.6 # orig: 0.0\n</weight>
  // Also need to catch id='X' (MG2.6.X) and id="X"
  std::regex rgx(R"(<weight.*id=[\"\'][^\d]*(\d+)[\"\'].*>([\s\S]*)</weight>)");
  for (auto it = lhe_info->headers_begin(); it != lhe_info->headers_end();
       ++it) {
    std::vector<std::string>::const_iterator iLt = it->begin();
    for (; iLt != it->end(); ++iLt) {
      std::string line = *iLt;
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
        continue;
      }
      if (line.find("</weightgroup") != std::string::npos) {
        lheWeightLabels_.push_back(line);
        lheWeightWasKept_.push_back(false);
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

          bool keep = false;
          for (auto const& x : split_label) {
            auto rule = getVarRule("lheweights:" + x);
            if (!rule.zeroed) {
              unsigned id = boost::lexical_cast<int>(rgx_match[1]);
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

  // edm::Handle<double> jets_rho_handle;
  // if (do_jets_rho_) {
  //   event.getByLabel(input_jets_rho_, jets_rho_handle);
  //   info_->set_jet_rho(*jets_rho_handle);
  // }
  // edm::Handle<double> lepton_rho_handle;
  // if (do_leptons_rho_) {
  //   event.getByLabel(input_leptons_rho_, lepton_rho_handle);
  //   info_->set_lepton_rho(*lepton_rho_handle);
  // }

  // for (unsigned i = 0; i < weights_.size(); ++i) {
  //   edm::Handle<double> weight;
  //   event.getByLabel(weights_[i].second, weight);
  //   double weights_result = (*weight);
  //   info_->set_weight(weights_[i].first, weights_result);
  // }


  // for (unsigned i = 0; i < gen_weights_.size(); ++i) {
  //   edm::Handle<GenFilterInfo> weight;
  //   event.getByLabel(gen_weights_[i].second, weight);
  //   double weights_result = weight->filterEfficiency();
  //   info_->set_weight(gen_weights_[i].first, weights_result);
  // }

  // edm::Handle<edm::View<reco::Vertex> > vtxs_handle;
  // if (do_vertex_count_) {
  //   event.getByLabel(input_vertices_, vtxs_handle);
  //   info_->set_good_vertices(vtxs_handle->size());
  // }

  if (includeLHEWeights_) {
    edm::Handle<LHEEventProduct> lhe_handle;
    event.getByToken(lheToken_, lhe_handle);
    //double nominalLHEWeight = lhe_handle->hepeup().XWGTUP;
    double nominalLHEWeight = lhe_handle->weights()[0].wgt ;
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

  // info->Print(1);

  // if(do_lhe_weights_ || do_ht_){
  //   event.getByLabel(lhe_collection_, lhe_handle);
  //  /* Accessing global evt weights directly from the LHEEventProduct
  //     disabled and replaced by taking them from the GenEventInfoProduct
  //     which should exist in all samples, independent of which generator was used
  //     if(lhe_handle->hepeup().XWGTUP>=0){
  //     info_->set_weight("wt_mc_sign",1);
  //     } else info_->set_weight("wt_mc_sign",-1);
  //   */
  //   if (do_ht_){
  //     std::vector<lhef::HEPEUP::FiveVector> lheParticles = lhe_handle->hepeup().PUP;
  //     double lheHt = 0.;
  //     unsigned nOutgoingPartons = 0;
  //     std::vector<ROOT::Math::PxPyPzEVector> zll_cands;
  //     for(size_t idxPart = 0; idxPart < lheParticles.size();++idxPart){
  //      unsigned absPdgId = TMath::Abs(lhe_handle->hepeup().IDUP[idxPart]);
  //      unsigned status = lhe_handle->hepeup().ISTUP[idxPart];
  //      if(status==1 &&((absPdgId >=1 &&absPdgId<=6) || absPdgId == 21)){
  //        lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxPart][0],2) + TMath::Power(lheParticles[idxPart][1],2));
  //        nOutgoingPartons++;
  //       }
  //      if(status == 1 && (absPdgId ==11 || absPdgId == 13 || absPdgId ==15)){
  //        zll_cands.push_back(ROOT::Math::PxPyPzEVector(lheParticles[idxPart][0],lheParticles[idxPart][1],lheParticles[idxPart][2],lheParticles[idxPart][3]));
  //      }
  //     }
  //     if(zll_cands.size() == 2){
  //       info_->set_gen_mll((zll_cands[0]+zll_cands[1]).M());
  //     }
  //     info_->set_gen_ht(lheHt);
  //     info_->set_n_outgoing_partons(nOutgoingPartons);
  //   }
  //   if (do_lhe_weights_) {
  //     double nominal_wt = lhe_handle->hepeup().XWGTUP;
  //     for (unsigned i = 0; i < lhe_handle->weights().size(); ++i) {
  //       info_->set_weight(lhe_handle->weights()[i].id,
  //                         lhe_handle->weights()[i].wgt / nominal_wt, false);
  //     }
  //   }
  // }

  // for (unsigned i = 0; i < filters_.size(); ++i) {
  //   edm::Handle<bool> filter;
  //   event.getByLabel(filters_[i].second, filter);
  //   bool filter_result = (*filter);
  //   if (invert_filter_logic_.find(filters_[i].first) !=
  //       invert_filter_logic_.end())
  //     filter_result = !filter_result;
  //   info_->set_filter_result(filters_[i].first, filter_result);
  //   observed_filters_[filters_[i].first] = CityHash64(filters_[i].first);
  // }

  // if (do_filtersfromtrig_) {
  //   edm::Handle<edm::TriggerResults> triggerResults;
  //   event.getByLabel(filtersfromtrig_input_, triggerResults);
  //   if (!triggerResults.isValid()) {
  //     throw cms::Exception("TriggerNotValid")
  //         << "Trigger Results is not valid\n";
  //   }

  //   const edm::TriggerNames& triggerNames = event.triggerNames(*triggerResults);
  //   for (unsigned iTrigger = 0; iTrigger < triggerResults->size(); iTrigger++) {
  //     std::string trigName = triggerNames.triggerName(iTrigger);
  //     bool pass = false;
  //     if (filtersfromtrig_.size() == 1 && filtersfromtrig_[0] == "*") {
  //       // copy all filters
  //       pass = true;
  //     } else {
  //       // copy only filters given in input
  //       for (unsigned iFilter = 0; iFilter < filtersfromtrig_.size();
  //            iFilter++) {
  //         std::string filtername = filtersfromtrig_[iFilter];
  //         if (filtername.find(trigName) != std::string::npos) {
  //           pass = true;
  //           break;
  //         }
  //       }
  //       // if (!pass) {
  //       // throw cms::Exception("TriggerNotFound")<<trigName<<" was not found in
  //       // trigger results\n";
  //       //}
  //     }
  //     if (!pass) continue;
  //     if (!triggerResults->wasrun(iTrigger)) {
  //       throw cms::Exception("TriggerNotRun") << trigName << " was not run\n";
  //     }
  //     bool filter_result = triggerResults->accept(iTrigger);
  //     if (invert_filter_logic_.find(trigName) != invert_filter_logic_.end())
  //       filter_result = !filter_result;

  //     // std::cout << " Filter " << trigName << " " <<
  //     // triggerResults->accept(iTrigger) << " " << filter_result << std::endl;

  //     info_->set_filter_result(trigName, filter_result);
  //     observed_filters_[trigName] = CityHash64(trigName);
  //   }  // loop on triggers
  // }

  // edm::Handle<reco::BeamHaloSummary> beam_halo_handle;
  // if (do_csc_filter_) {
  //   event.getByLabel(input_csc_filter_, beam_halo_handle);
  //   info_->set_filter_result("CSCTightHaloFilter",
  //                            beam_halo_handle->CSCTightHaloId());
  //   observed_filters_["CSCTightHaloFilter"] = CityHash64("CSCTightHaloFilter");
  //  }
}

void AcornEventInfoProducer::endStream() {
  // if (!observed_filters_.empty()) {
  //   std::cout << std::string(78, '-') << "\n";
  //   std::cout << boost::format("%-56s  %20s\n") %
  //                    std::string("EventInfo Filters") %
  //                    std::string("Hash Summmary");
  //   std::map<std::string, std::size_t>::const_iterator iter;
  //   for (iter = observed_filters_.begin(); iter != observed_filters_.end();
  //        ++iter) {
  //     std::cout << boost::format("%-56s| %020i\n") % iter->first % iter->second;
  //   }
  // }
  if (lheWeightLabels_.size()) {
    std::cout << "LHE event weights\n";
    std::cout << std::string(78, '-') << "\n";
    for (unsigned l = 0; l < lheWeightLabels_.size(); ++l) {
      if (lheWeightWasKept_[l]) std::cout << "(*) ";
      std::cout << lheWeightLabels_[l];
    }
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornEventInfoProducer);
