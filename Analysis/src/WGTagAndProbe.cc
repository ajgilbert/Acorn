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
#include "Acorn/Analysis/interface/WGTagAndProbe.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/Analysis/interface/WGAnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/Reduction.h"

namespace ac {

WGTagAndProbe::WGTagAndProbe(std::string const& name)
    : ModuleBase(name),
      fs_(nullptr),
      year_(2016),
      is_data_(true),
      do_photons_(false) {}

WGTagAndProbe::~WGTagAndProbe() { ; }

int WGTagAndProbe::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("WGTagAndProbe", "WGTagAndProbe");
    tree_->Branch("run", &run_);
    tree_->Branch("n_vtx", &n_vtx_);
    tree_->Branch("t_pt", &t_pt_);
    tree_->Branch("t_eta", &t_eta_);
    tree_->Branch("t_phi", &t_phi_);
    tree_->Branch("t_q", &t_phi_);
    tree_->Branch("t_id", &t_id_);
    tree_->Branch("t_rand", &t_rand_);

    tree_->Branch("p_pt", &p_pt_);
    tree_->Branch("p_eta", &p_eta_);
    tree_->Branch("p_phi", &p_phi_);
    tree_->Branch("p_q", &p_phi_);
    tree_->Branch("p_id", &p_id_);

    tree_->Branch("m_ll", &m_ll_);

    tree_->Branch("t_trg", &t_trg_);
    tree_->Branch("p_trg", &p_trg_);
    tree_->Branch("wt_def", &wt_def_);
  }

  if (is_data_) {
    filters_Ele27_ =
        LookupFilter({{272023, "hltEle27WPTightGsfTrackIsoFilter"}});

    filters_Ele32_ =
        LookupFilter({{281010, "hltEle32noerWPTightGsfTrackIsoFilter"},
                      {302023, "hltEle32WPTightGsfTrackIsoFilter"}});

    filters_Ele32_L1DoubleEG_ =
        LookupFilter({{295982, "hltEle32L1DoubleEGWPTightGsfTrackIsoFilter"}});

    filters_Ele32_L1DoubleEG_seed_ = // this is the wrong filter!
        LookupFilter({{295982, "hltL1sSingleAndDoubleEGor"}});

  } else {
    filters_Ele27_ =
        LookupFilter({{2016, "hltEle27WPTightGsfTrackIsoFilter"},
                      {2017, "hltEle27WPTightGsfTrackIsoFilter"},
                      {2018, "hltEle27WPTightGsfTrackIsoFilter"}});

    filters_Ele32_ =
        LookupFilter({{2016, "hltEle32noerWPTightGsfTrackIsoFilter"},
                      {2017, "hltEle32WPTightGsfTrackIsoFilter"},
                      {2018, "hltEle32WPTightGsfTrackIsoFilter"}});
  }

  TFile f(corrections_.c_str());
  ws_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  f.Close();
  fns_["pileup_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("pileup_ratio")->functor(ws_->argSet("pu_int")));
  return 0;
  }

  int WGTagAndProbe::Execute(TreeEvent* event) {

    SetDefaults();

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    run_ = info->run();
    n_vtx_ = info->numVertices();

    // auto muons = event->GetPtrVec<ac::Muon>("muons");
    auto electrons = event->GetPtrVec<ac::Electron>("electrons");
    auto photons = event->GetPtrVec<ac::Photon>("photons");

    for (Photon *p : photons) {
      PhotonIsoCorrector(p, info->userDoubles().at(0));
    }

    auto tags = ac::copy_keep_if(electrons, [](ac::Electron *e) {
      return e->pt() > 35. && fabs(e->scEta()) < 2.5 && e->isCutBasedMediumElectron() &&
                   (fabs(e->scEta()) < 1.4442 || fabs(e->scEta()) > 1.566) && ElectronIPCuts(e);
    });

    if (!is_data_) {
      wt_def_ = info->totalWeight();
      auto const& pu_info = event->GetPtrVec<PileupInfo>("pileupInfo");
      for (PileupInfo const* pu : pu_info) {
        if (pu->bunchCrossing() == 0) {
          wt_def_ *= RooFunc(fns_["pileup_ratio"], {pu->trueNumInteractions()});
          break;
        }
      }
    }


    if (!do_photons_) {
      auto probes = ac::copy_keep_if(electrons, [](ac::Electron *e) {
        return e->pt() > 35. && fabs(e->scEta()) < 2.5 &&
                     (fabs(e->scEta()) < 1.4442 || fabs(e->scEta()) > 1.566);
      });

      auto electron_pairs = ac::MakePairs(tags, probes);
      ac::keep_if(electron_pairs, [](std::pair<ac::Electron*, ac::Electron*> const& p) {
        return (p.first->charge() * p.second->charge() == -1) && ac::DeltaR(p.first, p.second) > 0.4;
      });


      bool req_pos = rng.Uniform() > 0.5;

      for (auto const& pair : electron_pairs) {
        ac::Electron const* ele_t = pair.first;
        ac::Electron const* ele_p = pair.second;

        t_pt_ = ele_t->pt();
        t_eta_ = ele_t->scEta();
        t_phi_ = ele_t->phi();
        t_q_ = ele_t->charge();
        t_id_ = ele_t->isCutBasedMediumElectron() && ElectronIPCuts(ele_t);
        t_rand_ = (ele_t->charge() == +1) == req_pos;

        p_pt_ = ele_p->pt();
        p_eta_ = ele_p->scEta();
        p_phi_ = ele_p->phi();
        p_q_ = ele_p->charge();
        p_id_ = ele_p->isCutBasedMediumElectron() && ElectronIPCuts(ele_p);

        m_ll_ = (ele_t->vector() + ele_p->vector()).M();

        t_trg_ = PassesTrigger(ele_t, event);
        p_trg_ = PassesTrigger(ele_p, event);

        tree_->Fill();
      }
    } else {
      auto probes = ac::copy_keep_if(photons, [&](ac::Photon const* p) {
      return p->pt() > 30. && fabs(p->scEta()) < 2.5 &&
             (fabs(p->scEta()) < 1.4442 || fabs(p->scEta()) > 1.566) && PhotonIDIso(p, year_, 1, true, true);
      });

      auto electron_pairs = ac::MakePairs(tags, probes);
      ac::keep_if(electron_pairs, [](std::pair<ac::Electron*, ac::Photon*> const& p) {
        return ac::DeltaR(p.first, p.second) > 0.4;
      });


      bool req_pos = rng.Uniform() > 0.5;

      for (auto const& pair : electron_pairs) {
        ac::Electron const* ele_t = pair.first;
        ac::Photon const* ele_p = pair.second;

        t_pt_ = ele_t->pt();
        t_eta_ = ele_t->scEta();
        t_phi_ = ele_t->phi();
        t_q_ = ele_t->charge();
        t_id_ = ele_t->isCutBasedMediumElectron() && ElectronIPCuts(ele_t);
        t_rand_ = (ele_t->charge() == +1) == req_pos;

        p_pt_ = ele_p->pt();
        p_eta_ = ele_p->scEta();
        p_phi_ = ele_p->phi();
        p_q_ = ele_p->charge();
        p_id_ = (ele_p->worstChargedIsolation() - ele_p->chargedIso()) < std::min(0.05 * ele_p->pt(), 6.0);

        m_ll_ = (ele_t->vector() + ele_p->vector()).M();

        t_trg_ = PassesTrigger(ele_t, event);
        p_trg_ = true;

        tree_->Fill();
      }
    }


    return 0;
  }

  void WGTagAndProbe::SetDefaults() {
    run_ = 0;
    n_vtx_ = 0;
    t_pt_ = 0.;
    t_eta_ = 0.;
    t_phi_ = 0.;
    t_q_ = 0;
    t_id_ = false;
    p_pt_ = 0.;
    p_eta_ = 0.;
    p_phi_ = 0.;
    p_id_ = false;
    p_q_ = 0;
    m_ll_ = 0.;
    t_trg_ = false;
    p_trg_ = false;
    wt_def_ = 1.;
  }

  int WGTagAndProbe::PostAnalysis() {
    return 0;
  }

  void WGTagAndProbe::PrintInfo() {}

  bool WGTagAndProbe::PassesTrigger(ac::Electron const* e, ac::TreeEvent * event) const {
    unsigned trg_lookup = is_data_ ? run_ : year_;
    bool e0_trg = false;
    if (year_ == 2016) {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele27_WPTight_Gsf");
      e0_trg = IsFilterMatchedDR(e, trg_objs, filters_Ele27_.Lookup(trg_lookup), 0.3);
    } else if (year_ == 2017) {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf");
      e0_trg = IsFilterMatchedDR(e, trg_objs, filters_Ele32_.Lookup(trg_lookup), 0.3);
      if (is_data_ && run_ < 302026) {
        auto const& trg_objs_alt = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf_L1DoubleEG");
        bool l0_trg_alt = IsFilterMatchedDR(e, trg_objs_alt, filters_Ele32_L1DoubleEG_.Lookup(trg_lookup), 0.3) &&
                          IsFilterMatchedDR(e, trg_objs_alt, filters_Ele32_L1DoubleEG_seed_.Lookup(trg_lookup), 0.3);
        e0_trg = e0_trg || l0_trg_alt;
      }
    } else if (year_ == 2018) {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf");
      e0_trg = IsFilterMatchedDR(e, trg_objs, filters_Ele32_.Lookup(trg_lookup), 0.3);
    }
    return e0_trg;
  }

}
