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

void JESHists::MakeHists(std::vector<std::string> & srcs) {
  for (auto src : srcs) {
    jhists_.push_back(dir.make<TH2D>(src.c_str(), "", pt_edges_.size() - 1, pt_edges_.data(), eta_edges_.size() - 1, eta_edges_.data()));
  }
    j_nominal_ = dir.make<TH2D>("nominal", "", pt_edges_.size() - 1, pt_edges_.data(), eta_edges_.size() - 1, eta_edges_.data());
}

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


  for (auto src : jsrcs_) {
    auto pars = new JetCorrectorParameters(jfile_, src);
    jpars_.push_back(pars);
    juncs_.push_back(new JetCorrectionUncertainty(*pars));
  }

  auto & main = jhists_["main"];
  main.dir = fs_->mkdir("main");
  main.pt_edges_ = {40, 60, 80, 100, 150, 200};
  main.eta_edges_ = {-5.0, -3.0, -2.5, -1.3, 0.0, +1.3, +2.5, +3.0, +5.0};
  main.MakeHists(jsrcs_);

  auto & fine = jhists_["fine"];
  fine.dir = fs_->mkdir("fine");
  fine.pt_edges_ = {35, 40, 45, 50, 60, 70, 80, 100, 150, 200, 300};
  fine.eta_edges_ = {-5.0, -3.0, -2.5, -1.3, 0.0, +1.3, +2.5, +3.0, +5.0};
  fine.MakeHists(jsrcs_);

  auto & fineneg = jhists_["fineneg"];
  fineneg.dir = fs_->mkdir("fineneg");
  fineneg.pt_edges_ = {35, 40, 45, 50, 60, 70, 80, 100, 150, 200, 300};
  fineneg.eta_edges_ = {-5.0, -3.0, -2.5, -1.3, 0.0, +1.3, +2.5, +3.0, +5.0};
  fineneg.MakeHists(jsrcs_);

  auto & direct = jhists_["direct"];
  direct.dir = fs_->mkdir("direct");
  // direct.pt_edges_ = {30, 35, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300};
  // direct.eta_edges_ = {-5.0, -3.0, -2.5, -1.3, 0.0, +1.3, +2.5, +3.0, +5.0};
  direct.pt_edges_ = {40, 60, 80, 100, 150, 200};
  direct.eta_edges_ = {-5.0, -3.0, -2.5, -1.3, 0.0, +1.3, +2.5, +3.0, +5.0};
  // direct.eta_edges_ = {-5.0, -2.5, 0.0, +2.5, +5.0};
  direct.MakeHists(jsrcs_);

  for (int ix = 1; ix <= direct.j_nominal_->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= direct.j_nominal_->GetNbinsY(); ++iy) {
      direct.j_nominal_->SetBinContent(ix, iy, 1.0);
      for (unsigned ic = 0; ic < jsrcs_.size(); ++ic) {
        juncs_[ic]->setJetPt(direct.j_nominal_->GetXaxis()->GetBinCenter(ix));
        juncs_[ic]->setJetEta(direct.j_nominal_->GetYaxis()->GetBinCenter(iy));
        direct.jhists_[ic]->SetBinContent(ix, iy, 1. + juncs_[ic]->getUncertainty(true));
      }
    }
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

  auto & main = jhists_["main"];
  auto & fine = jhists_["fine"];
  auto & fineneg = jhists_["fineneg"];
  for (auto const& j : jets) {
    main.j_nominal_->Fill(j->pt(), j->eta());
    fine.j_nominal_->Fill(j->pt(), j->eta());
    fineneg.j_nominal_->Fill(j->pt(), j->eta());
    for (unsigned ic = 0; ic < jsrcs_.size(); ++ic) {
      juncs_[ic]->setJetPt(j->pt());
      juncs_[ic]->setJetEta(j->eta());
      double uncHi = juncs_[ic]->getUncertainty(true);

      main.jhists_[ic]->Fill(j->pt() * (1. + uncHi), j->eta());
      fine.jhists_[ic]->Fill(j->pt() * (1. + uncHi), j->eta());
      fineneg.jhists_[ic]->Fill(j->pt() * (1. - uncHi), j->eta());
    }
  }

  return 0;
}

int JESStudy::PostAnalysis() { return 0; }

void JESStudy::PrintInfo() {}

}
