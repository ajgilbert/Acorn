#ifndef Acorn_NTupler_AcornBaseProducer_h
#define Acorn_NTupler_AcornBaseProducer_h

#include <memory>
#include <regex>
#include <string>
#include <vector>
#include "Acorn/NTupler/interface/Reduction.h"
#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Rtypes.h"
#include "boost/algorithm/string.hpp"
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"

template <class T>
class AcornBaseProducer : public edm::stream::EDProducer<> {
 public:
  explicit AcornBaseProducer(const edm::ParameterSet &config)
      : branch_(config.getParameter<std::string>("branch")) {
    // Configure the list of keep/drop statements

    auto save_config = config.getParameter<std::vector<std::string>>("select");
    save_config_.push_back("keep .*");
    save_config_.insert(std::end(save_config_), std::begin(save_config), std::end(save_config));
    for (auto cfg : save_config_) {
      // remove leading or trailing whitespace
      boost::trim(cfg);

      // split line into words, allowing for any amount of whitespace
      std::vector<std::string> split_cfg;
      boost::split(split_cfg, cfg, boost::is_any_of(" \t"), boost::token_compress_on);

      // Throw this exception later if any problem parsing the config is encountered
      auto except = cms::Exception("ConfigError")
              << "The select statement: \"" << cfg << "\" could not be parsed\n";

      // Validation: must have at least two words, the first of which is either "keep" or "drop"
      if (split_cfg.size() < 2 || !(split_cfg.at(0) == "keep" || split_cfg.at(0) == "drop")) {
        throw except;
      }

      for (unsigned i = 1; i < split_cfg.size(); ++i) {
        VarRule newrule;
        newrule.zeroed = split_cfg.at(0) == "drop" ? true : false;
        std::vector<std::string> split_pattern;
        boost::split(split_pattern, split_cfg.at(i), boost::is_any_of("="));
        if (!(split_pattern.size() == 1 || split_pattern.size() == 2) ||
            split_pattern.at(0).size() == 0) {
          throw except;
        }
        newrule.rgx = std::regex(split_pattern.at(0));
        newrule.original = split_pattern.at(0);
        if (split_pattern.size() == 2) {
          newrule.truncate = true;
          newrule.round = boost::lexical_cast<int>(split_pattern.at(1));
        }
        var_rules_.push_back(newrule);
      }
    }

    output_ = new T();
  }

  virtual ~AcornBaseProducer() { delete output_; }

  virtual void beginStream(edm::StreamID id) { streamid_ = id.value(); }

  virtual void endStream() {
    // std::cout << "Producer for branch: " << branch_ << "\n";
    // for (auto const& it : resolved_rules_) {
    //   std::cout << it.first << ": " << "\t" << it.second.zeroed << "\t" << it.second.truncate << "\t" << it.second.round << "\n";
    // }
   }

  virtual void beginRun(edm::Run const &, edm::EventSetup const &) {
    if (!attachedBranch_) {
      // AcornEventProducer::getStreamTree(streamid_)->Branch(branch_.c_str(), &output_);
      AcornEventProducer::AddBranch(streamid_, branch_, output_);
    }
    attachedBranch_ = true;
  }


  struct VarRule {
    bool zeroed = false;
    std::regex rgx;
    std::string original = "";
    bool truncate = false;
    int round = 0;
  };

  VarRule getVarRule(std::string const& label) {
    VarRule resolved_rule;
    // Check the list of rules to see what we're doing with this var,
    // otherwise just use the default
    for (auto const &rule : var_rules_) {
      if (std::regex_match(label, rule.rgx)) {
        resolved_rule = rule;
      }
    }
    return resolved_rule;
  }

  template <typename D>
  D processVar(D const &arg, VarRule const& rule) {
    return rule.zeroed ? D() : (rule.truncate ? Reduce<D>(arg, rule.round) : arg);
  }

  template <typename D>
  D setVar(std::string const &label, D const &arg) {
    auto it = resolved_rules_.find(label);
    if (it == resolved_rules_.end()) {
      VarRule resolved_rule = getVarRule(label);
      it = resolved_rules_.emplace(label, resolved_rule).first;
    }
    return processVar(arg, it->second);
  }


  inline T *output() { return output_; }


 private:
  virtual void produce(edm::Event &, const edm::EventSetup &) = 0;

  std::vector<std::string> save_config_;      // the list of keep/save statements given by the user
  std::vector<VarRule> var_rules_;            // statements parsed into rules
  std::map<std::string, VarRule> resolved_rules_;  // resolved rules for specific variable names
  T *output_;  // pointer to the output object that will be added to the TTree
  std::string branch_;  // branch name
  unsigned streamid_;
  bool attachedBranch_ = false;
};

#endif
