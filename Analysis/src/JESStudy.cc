#include "Acorn/Analysis/interface/JESStudy.h"
#include <algorithm>
#include <map>
#include "TH2D.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/Analysis/interface/WGAnalysisTools.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/PFJet.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/Reduction.h"
#include "Math/Boost.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Rotation3D.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TMath.h"
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/range/algorithm/sort.hpp"

namespace ac {

JESStudy::JESStudy(std::string const& name)
    : ModuleBase(name), fs_(nullptr), year_(2016), is_data_(true) {}

JESStudy::~JESStudy() { ; }

int JESStudy::PreAnalysis() {
  jsrcs_ = {
    "AbsoluteStat",
    "AbsoluteScale",
    "AbsoluteMPFBias",
    "AbsoluteSample",
    "Fragmentation",
    "SinglePionECAL",
    "SinglePionHCAL",
    "FlavorQCD",
    "TimePtEta",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeBal",
    "RelativeSample",
    "RelativeFSR",
    "RelativeStatFSR",
    "RelativeStatEC",
    "RelativeStatHF",
    "PileUpDataMC",
    "PileUpPtRef",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF"
  };

  std::string jfile_ = "wgamma/inputs/Autumn18_V8_MC_UncertaintySources_AK4PFchs.txt";

  std::vector<double> pt_edges = {40, 60, 80, 100, 150, 200};
  std::vector<double> eta_edges = {-5.0, -3.0, -2.5, -1.3, 0.0, +1.3, +2.5, +3.0, +5.0};
  for (auto src : jsrcs_) {
    auto pars = new JetCorrectorParameters(jfile_, src);
    jpars_.push_back(pars);
    juncs_.push_back(new JetCorrectionUncertainty(*pars));
    jhists_.push_back(fs_->make<TH2D>(src.c_str(), "", pt_edges.size() - 1, pt_edges.data(), eta_edges.size() - 1, eta_edges.data()));
  }

  if (fs_) {
    // tree_ = fs_->make<TTree>("JESStudy", "JESStudy");
    j_nominal = fs_->make<TH2D>("nominal", "", pt_edges.size() - 1, pt_edges.data(), eta_edges.size() - 1, eta_edges.data());
  }
  return 0;
}

int JESStudy::Execute(TreeEvent* event) {
  // auto const* info = event->GetPtr<EventInfo>("eventInfo");

  auto all_muons = event->GetPtrVec<ac::Muon>("muons");
  auto all_elecs = event->GetPtrVec<ac::Electron>("electrons");

  auto muons = ac::copy_keep_if(all_muons, [](ac::Muon const* m) {
    return m->pt() > 20. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && MuonPFIso(m) < 0.15 &&
           fabs(m->dxy()) < 0.05 && fabs(m->dz()) < 0.2;
  });

  auto elecs = ac::copy_keep_if(all_elecs, [](ac::Electron const* e) {
    return e->pt() > 20. && fabs(e->scEta()) < 2.5 && e->isCutBasedMediumElectron() &&
           ElectronIsoFall17V2(e, 2) && (fabs(e->scEta()) < 1.4442 || fabs(e->scEta()) > 1.566) &&
           ElectronIPCuts(e);
  });
  boost::range::sort(muons, DescendingPt);
  boost::range::sort(elecs, DescendingPt);

  bool is_zmm = (muons.size() == 2 && muons[0]->charge() != muons[1]->charge() && DeltaR(muons[0], muons[1]) > 0.3);
  bool is_zee = (elecs.size() == 2 && elecs[0]->charge() != elecs[1]->charge() && DeltaR(elecs[0], elecs[1]) > 0.3);

  if (!(is_zmm || is_zee)) return 0;

  auto all_jets = event->GetPtrVec<ac::PFJet>("pfJets");
  ac::Candidate *l0 = is_zmm ? (ac::Candidate*)muons[0] : (ac::Candidate*)elecs[0];
  ac::Candidate *l1 = is_zmm ? (ac::Candidate*)muons[1] : (ac::Candidate*)elecs[1];

  auto jets = ac::copy_keep_if(all_jets, [&](ac::PFJet const* j) {
    return DeltaR(j, l0) > 0.5 && DeltaR(j, l1) > 0.5 && j->passesJetID();
  });

  for (auto const& j : jets) {
    j_nominal->Fill(j->pt(), j->eta());
    for (unsigned ic = 0; ic < jsrcs_.size(); ++ic) {
      juncs_[ic]->setJetPt(j->pt());
      juncs_[ic]->setJetEta(j->eta());
      jhists_[ic]->Fill(j->pt() * (1. + juncs_[ic]->getUncertainty(true)), j->eta());
    }
  }

  return 0;
}

int JESStudy::PostAnalysis() { return 0; }

void JESStudy::PrintInfo() {}

}
