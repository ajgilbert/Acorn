#include <algorithm>
#include <map>
#include "TMath.h"
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/range/algorithm/sort.hpp"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "Acorn/Analysis/interface/WGDataAnalysis.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"

namespace ac {

WGDataAnalysis::WGDataAnalysis(std::string const& name)
    : ModuleBase(name), fs_(nullptr), year_(2016), is_data_(true) {}

WGDataAnalysis::~WGDataAnalysis() { ; }

int WGDataAnalysis::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("WGDataAnalysis", "WGDataAnalysis");
    tree_->Branch("pt_m", &pt_m_);
    tree_->Branch("pt_p", &pt_p_);
    tree_->Branch("eta_m", &eta_m_);
    tree_->Branch("eta_p", &eta_p_);
    tree_->Branch("pt_met", &pt_met_);
    tree_->Branch("mt", &mt_);
    tree_->Branch("trg_m", &trg_m_);
    tree_->Branch("wt_pu", &wt_pu_);
    tree_->Branch("wt_m", &wt_m_);
    tree_->Branch("wt_trg_m", &wt_trg_m_);
    tree_->Print();
  }

  if (is_data_) {
    filters_IsoMu24_ =
        LookupFilter({{272023, "hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"},
                      {295982, "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"}});

    filters_IsoTkMu24_ =
        LookupFilter({{272023, "hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"},
                      {295982, "hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p07"}});

    filters_IsoMu27_ =
        LookupFilter({{272023, "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"},
                      {295982, "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"}});
  } else {
    filters_IsoMu24_ =
        LookupFilter({{2016, "hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"},
                      {2017, "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"},
                      {2018, "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"}});

    filters_IsoTkMu24_ =
        LookupFilter({{2016, "hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"}});

    filters_IsoMu27_ =
        LookupFilter({{2016, "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"},
                      {2017, "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"},
                      {2018, "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"}});
  }

  TFile f(corrections_.c_str());
  ws_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  f.Close();
  fns_["pileup_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("pileup_ratio")->functor(ws_->argSet("pu_int")));
  fns_["m_idisotrk_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_idisotrk_ratio")->functor(ws_->argSet("m_pt,m_eta")));
  fns_["m_trg_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_trg24_ratio")->functor(ws_->argSet("m_pt,m_eta")));

  return 0;
  }

  int WGDataAnalysis::Execute(TreeEvent* event) {

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    auto muons = event->GetPtrVec<ac::Muon>("muons");
    auto photons = event->GetPtrVec<ac::Photon>("photons");
    auto met = event->GetPtrVec<ac::Met>("pfType1Met");

    // Apply pT and ID cuts
    // The medium ID already includes iso?
    ac::keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && MuonPFIso(m) < 0.15;
    });

    ac::keep_if(photons, [](ac::Photon const* p) {
      return p->pt() > 30. && fabs(p->eta()) < 2.5 && p->isMediumIdPhoton() && !p->hasPixelSeed();
    });

    boost::range::sort(muons, DescendingPt);
    boost::range::sort(photons, DescendingPt);

    if (!(muons.size() == 1 && photons.size() == 1)) {
      return 1;
    }

    double dr_m_p = ac::DeltaR(muons[0], photons[0]);

    if (!(dr_m_p > 0.5)) {
      return 1;
    }

    if (muons.size() == 1 && photons.size() == 1) {

      pt_m_ = muons[0]->pt();
      eta_m_ = muons[0]->eta();

      pt_p_ = photons[0]->pt();
      eta_p_ = photons[0]->eta();

      pt_met_ = met[0]->pt();

      mt_ = ac::MT(muons[0], met[0]);

      unsigned trg_lookup = is_data_ ? info->run() : year_;
      if (year_ == 2016) {
        auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
        auto const& trg_objs_tk = event->GetPtrVec<TriggerObject>("triggerObjects_IsoTkMu24");
        trg_m_ =
            IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3) ||
            IsFilterMatchedDR(muons[0], trg_objs_tk, filters_IsoTkMu24_.Lookup(trg_lookup), 0.3);
      } else {
        auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu27");
        trg_m_ = IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu27_.Lookup(trg_lookup), 0.3);
      }

      wt_pu_ = 1.;
      wt_m_ = 1.;
      wt_trg_m_ = 1.;

      if (!is_data_) {
        auto const& pu_info = event->GetPtrVec<PileupInfo>("pileupInfo");
        for (PileupInfo const* pu : pu_info) {
          if (pu->bunchCrossing() == 0) {
            wt_pu_ = RooFunc(fns_["pileup_ratio"], {pu->trueNumInteractions()});
            break;
          }
        }
        wt_m_ = RooFunc(fns_["m_idisotrk_ratio"], {pt_m_, eta_m_});
        wt_trg_m_ = RooFunc(fns_["m_trg_ratio"], {pt_m_, eta_m_});
      }
      tree_->Fill();
    }


    return 0;
  }
  int WGDataAnalysis::PostAnalysis() {
    return 0;
  }

  void WGDataAnalysis::PrintInfo() {}



}
