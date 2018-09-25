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
  for (auto const& item: cfg.items()) {
    if (item.key() == "binned") {
      binned_ = item.value()[0];
    } else {
      json const& val = item.value();
      SampleInfo info;
      info.name = item.key();
      info.target = item.value()["target"];
      info.has_min = false;
      info.min = 0.;
      if (val["min"].size()) {
        info.has_min = true;
        info.min = val["min"][0];
      }
      info.has_max = false;
      info.max = 0.;
      if (val["max"].size()) {
        info.has_max = true;
        info.max = val["max"][0];
      }
      info.events = val["events"];
      info.xsec = val["xsec"];
      std::cout << info.name << "\t" << info.target << "\t" << info.has_min << "\t" << info.has_max
                << "\t" << info.min << "\t" << info.max << "\t" << info.events << "\t" << info.xsec
                << "\n";
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
  if (binned_ == "photon_pt") {
    var1_ = 0.0;
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
    var1_ = photon->pt();
  }

  // Then this part is generic
  double eff_lumi = 0.;
  double target_lumi = 0.;
  for (unsigned i = 0; i < samples_.size(); ++i) {
    bool in_sample = (((!samples_[i].has_min) || var1_ >= samples_[i].min) &&
                      ((!samples_[i].has_max) || var1_ < samples_[i].max));
    if (!in_sample) {
      continue;
    }
    double sample_lumi = (samples_[i].events / samples_[i].xsec);
    eff_lumi += sample_lumi;
    if (samples_[i].target) {
      target_lumi += sample_lumi;
    }
  }
  double weight =  target_lumi / eff_lumi;

  info->setWeight("stitching", weight);
  return 0;
}
int SampleStitching::PostAnalysis() { return 0; }

void SampleStitching::PrintInfo() {}
}  // namespace ac
