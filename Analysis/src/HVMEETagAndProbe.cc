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
#include "Acorn/Analysis/interface/HVMEETagAndProbe.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/Reduction.h"

namespace ac {

HVMEETagAndProbe::HVMEETagAndProbe(std::string const& name)
    : ModuleBase(name),
      fs_(nullptr),
      year_(2016),
      is_data_(true) {}

HVMEETagAndProbe::~HVMEETagAndProbe() { ; }

int HVMEETagAndProbe::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("HVMEETagAndProbe", "HVMEETagAndProbe");
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

    filters_Ele35_ =
        LookupFilter({{294927, "hltEle35noerWPTightGsfTrackIsoFilter"}});

    filters_Ele32_ =
        LookupFilter({{281010, "hltEle32noerWPTightGsfTrackIsoFilter"},
                      {302023, "hltEle32WPTightGsfTrackIsoFilter"}});


  } else {

    filters_Ele35_ =
        LookupFilter({{2017, "hltEle35noerWPTightGsfTrackIsoFilter"}});

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

  int HVMEETagAndProbe::Execute(TreeEvent* event) {

    SetDefaults();

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    run_ = info->run();
    n_vtx_ = info->numVertices();

    auto electrons = event->GetPtrVec<ac::Electron>("electrons");

    auto tags = ac::copy_keep_if(electrons, [](ac::Electron *e) {
      return e->pt() > 15. && fabs(e->eta()) < 2.1 && e->isMVAwp80Electron();
    });

    auto probes = ac::copy_keep_if(electrons, [](ac::Electron *e) {
      return e->pt() > 15. && fabs(e->eta()) < 2.1;
    });

    auto ele_pairs = ac::MakePairs(tags, probes);
    ac::keep_if(ele_pairs, [](std::pair<ac::Electron*, ac::Electron*> const& p) {
      return (p.first->charge() * p.second->charge() == -1) && ac::DeltaR(p.first, p.second) > 0.4;
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

    bool req_pos = rng.Uniform() > 0.5;

    for (auto const& pair : ele_pairs) {
      ac::Electron const* e_t = pair.first;
      ac::Electron const* e_p = pair.second;

      t_pt_ = e_t->pt();
      t_eta_ = e_t->eta();
      t_phi_ = e_t->phi();
      t_q_ = e_t->charge();
      t_id_ = e_t->isMVAwp80Electron();
      t_rand_ = (e_t->charge() == +1) == req_pos;

      p_pt_ = e_p->pt();
      p_eta_ = e_p->eta();
      p_phi_ = e_p->phi();
      p_q_ = e_p->charge();
      p_id_ = e_p->isMVAwp80Electron();

      m_ll_ = (e_t->vector() + e_p->vector()).M();

      t_trg_ = PassesTrigger(e_t, event);
      p_trg_ = PassesTrigger(e_p, event);

      tree_->Fill();
    }

    return 0;
  }

  void HVMEETagAndProbe::SetDefaults() {
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

  int HVMEETagAndProbe::PostAnalysis() {
    return 0;
  }

  void HVMEETagAndProbe::PrintInfo() {}

  bool HVMEETagAndProbe::PassesTrigger(ac::Electron const* e, ac::TreeEvent * event) const {
    unsigned trg_lookup = is_data_ ? run_ : year_;
    bool e_trg=false;
    if (year_ == 2016) {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele27_WPTight_Gsf");
      e_trg =
          IsFilterMatchedDR(e, trg_objs, filters_Ele27_.Lookup(trg_lookup), 0.3);
    } else if (year_ == 2017){
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele35_WPTight_Gsf");
      e_trg = IsFilterMatchedDR(e, trg_objs, filters_Ele35_.Lookup(trg_lookup), 0.3);
    } else {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf");
      e_trg = IsFilterMatchedDR(e, trg_objs, filters_Ele32_.Lookup(trg_lookup), 0.3);
    }
    return e_trg; 
  }

}
