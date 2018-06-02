#include <algorithm>
#include <map>
#include "TMath.h"
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"
#include "RooRealVar.h"
#include "Acorn/Analysis/interface/DiMuonAnalysis.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/EventInfo.h"

namespace ac {

  DiMuonAnalysis::DiMuonAnalysis(std::string const& name)
      : ModuleBase(name), fs_(nullptr){}

  DiMuonAnalysis::~DiMuonAnalysis() { ; }

  int DiMuonAnalysis::PreAnalysis() {
    if (fs_) {
      tree_ = fs_->make<TTree>("DiMuonAnalysis", "DiMuonAnalysis");
      tree_->Branch("pt_1", &pt_1_);
      tree_->Branch("pt_2", &pt_2_);
      tree_->Branch("m_ll", &m_ll_);
    }
    return 0;
  }

  int DiMuonAnalysis::Execute(TreeEvent* event) {

    std::vector<ac::Muon *> muons = event->GetPtrVec<ac::Muon>("muons");
    ac::keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon();
    });

    std::sort(muons.begin(), muons.end(), [](Muon const* m1, Muon const* m2) {
      return m1->pt() > m2->pt();
    });

    if (muons.size() == 2 && (muons[0]->charge() * muons[1]->charge()) == -1) {
      pt_1_ = muons[0]->pt();
      pt_2_ = muons[1]->pt();
      m_ll_ = (muons[0]->vector() + muons[1]->vector()).M();
      tree_->Fill();
    }

    return 0;
  }
  int DiMuonAnalysis::PostAnalysis() {
    return 0;
  }

  void DiMuonAnalysis::PrintInfo() {}



}
