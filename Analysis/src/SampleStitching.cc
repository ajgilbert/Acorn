#include "Acorn/Analysis/interface/SampleStitching.h"
#include <algorithm>
#include <map>
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/json.hpp"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "Math/Boost.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Rotation3D.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"

namespace ac {

SampleStitching::SampleStitching(std::string const& name, nlohmann::json const& cfg)
    : ModuleBase(name), fs_(nullptr) {
  using json = nlohmann::json;
  binned_ = cfg["binned"].get<std::vector<std::string>>();
  unsigned N = binned_.size();
  vars_.resize(N);
  for (auto const& item: cfg.items()) {
    if (item.key() == "binned") {
      continue;
    } else {
      json const& val = item.value();
      SampleInfo info;
      info.has_min.resize(N);
      info.has_max.resize(N);
      info.min.resize(N);
      info.max.resize(N);
      info.name = item.key();
      info.target = val["target"];

      for (unsigned i = 0; i < binned_.size(); ++i) {
        info.has_min[i] = false;
        info.has_max[i] = false;
        if (val["min"].size() >= N && !val["min"][i].is_null()) {
          info.has_min[i] = true;
          info.min[i] = val["min"][i];
        }
        if (val["max"].size() >= N && !val["max"][i].is_null()) {
          info.has_max[i] = true;
          info.max[i] = val["max"][i];
        }
      }

      info.events = val["events"];
      info.xsec = val["xsec"];
      std::cout << info.name << "\t" << info.target;
      for (unsigned i = 0; i < binned_.size(); ++i) {
        std::cout << "\t" << info.has_min[i] << "\t" << info.has_max[i] << "\t" << info.min[i]
                  << "\t" << info.max[i];
      }
      std::cout << "\t" << info.events << "\t" << info.xsec << "\n";
      samples_.push_back(info);
    }
  }
}

SampleStitching::~SampleStitching() { ; }

int SampleStitching::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("SampleStitching", "SampleStitching");
    // tree_->Branch("wt", &wt_);
    // tree_->Branch("n_jets", &n_jets_);
  }
  return 0;
}

int SampleStitching::Execute(TreeEvent* event) {
  auto lhe_parts = event->GetPtrVec<ac::GenParticle>("lheParticles");
  auto info = event->GetPtr<ac::EventInfo>("eventInfo");

  // Here we extract the variables for the different kinds of stitching we know about
  for (unsigned i = 0; i < binned_.size(); ++i) {
    if (binned_[i] == "photon_pt") {
      vars_[i] = 0.0;
      ac::GenParticle* photon = nullptr;
      for (auto const& p : lhe_parts) {
        if (p->pdgId() == 22) {
          photon = p;
          break;
        }
      }
      if (!photon) {
        throw std::runtime_error(
            "SampleStitching: binning is \"photon_pt\" but no LHE photon was found");
      }
      vars_[i] = photon->pt();
    } else if (binned_[i] == "njet") {
      vars_[i] = 0.5;
      for (auto const& p : lhe_parts) {
        unsigned apdgid = std::abs(p->pdgId());
        if (p->status() == 1 && ((apdgid >= 1 && apdgid <= 6) || apdgid == 21)) {
          vars_[i] += 1.0;
        }
      }
    } else if (binned_[i] == "ht") {
      vars_[i] = 0.0;
      for (auto const& p : lhe_parts) {
        unsigned apdgid = std::abs(p->pdgId());
        if (p->status() == 1 && ((apdgid >= 1 && apdgid <= 6) || apdgid == 21)) {
          vars_[i] += p->pt();
        }
      }
    } else {
      throw std::runtime_error(
          "SampleStitching: No implementation exists for binning variable \"" + binned_[i] + "\"");
    }
  }

  // Then this part is generic
  double eff_lumi = 0.;
  double target_lumi = 0.;

  // std::cout << "vars: ";
  // for (unsigned j = 0; j < binned_.size(); ++j) {
  //   std::cout << vars_[j] << "\t";
  // }
  // std::cout << "\n";

  for (unsigned i = 0; i < samples_.size(); ++i) {
    bool in_sample = true;
    for (unsigned j = 0; j < binned_.size(); ++j) {
      in_sample = in_sample && (((!samples_[i].has_min[j]) || vars_[j] >= samples_[i].min[j]) &&
                                ((!samples_[i].has_max[j]) || vars_[j] < samples_[i].max[j]));
    }
    if (!in_sample) {
      continue;
    }
    // std::cout << " - in sample " << samples_[i].name << "\n";
    double sample_lumi = (samples_[i].events / samples_[i].xsec);
    eff_lumi += sample_lumi;
    if (samples_[i].target) {
      target_lumi += sample_lumi;
    }
  }
  // std::cout << " => " << target_lumi << "/" << eff_lumi << "\n";
  double weight =  target_lumi / eff_lumi;

  info->setWeight("stitching", weight);
  return 0;
}
int SampleStitching::PostAnalysis() { return 0; }

void SampleStitching::PrintInfo() {}
}  // namespace ac
