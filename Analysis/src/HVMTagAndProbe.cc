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
#include "Acorn/Analysis/interface/HVMTagAndProbe.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/Reduction.h"

namespace ac {

HVMTagAndProbe::HVMTagAndProbe(std::string const& name)
    : ModuleBase(name),
      fs_(nullptr),
      year_(2016),
      is_data_(true) {}

HVMTagAndProbe::~HVMTagAndProbe() { ; }

int HVMTagAndProbe::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("HVMTagAndProbe", "HVMTagAndProbe");
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
  return 0;
  }

  int HVMTagAndProbe::Execute(TreeEvent* event) {

    SetDefaults();

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    run_ = info->run();
    n_vtx_ = info->numVertices();

    auto muons = event->GetPtrVec<ac::Muon>("muons");
    std::vector<ac::Track *> tracksforiso = event->GetPtrVec<ac::Track>("TracksForIso");
    //auto electrons = event->GetPtrVec<ac::Electron>("electrons");

    auto tags = ac::copy_keep_if(muons, [](ac::Muon *m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && MuonPFIso(m)<0.15;
    });

    auto probes = ac::copy_keep_if(muons, [](ac::Muon *m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4;
    });

    auto muon_pairs = ac::MakePairs(tags, probes);
    ac::keep_if(muon_pairs, [](std::pair<ac::Muon*, ac::Muon*> const& p) {
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

    for (auto const& pair : muon_pairs) {
      ac::Muon const* mu_t = pair.first;
      ac::Muon const* mu_p = pair.second;

      t_pt_ = mu_t->pt();
      t_eta_ = mu_t->eta();
      t_phi_ = mu_t->phi();
      t_q_ = mu_t->charge();
      t_id_ = mu_t->isMediumMuon() && MuonPFIso(mu_t)<0.15;
      t_rand_ = (mu_t->charge() == +1) == req_pos;

      p_pt_ = mu_p->pt();
      p_eta_ = mu_p->eta();
      p_phi_ = mu_p->phi();
      p_q_ = mu_p->charge();
      p_id_ = mu_p->isMediumMuon() && MuonPFIso(mu_p)<0.15;

      m_ll_ = (mu_t->vector() + mu_p->vector()).M();

      t_trg_ = PassesTrigger(mu_t, event);
      p_trg_ = PassesTrigger(mu_p, event);

      p_trk_iso_=0;
       std::cout<<"for this pair"<<std::endl;
       for (unsigned i=0; i <tracksforiso.size(); i++){
         if(DeltaRTrack(tracksforiso.at(i),mu_p)<0.3) p_trk_iso_+=tracksforiso.at(i)->pt();
         std::cout<<DeltaRTrack(tracksforiso.at(i),mu_p)<<std::endl;
       }

      tree_->Fill();
    }

    return 0;
  }

  void HVMTagAndProbe::SetDefaults() {
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

  int HVMTagAndProbe::PostAnalysis() {
    return 0;
  }

  void HVMTagAndProbe::PrintInfo() {}

  bool HVMTagAndProbe::PassesTrigger(ac::Muon const* m, ac::TreeEvent * event) const {
    unsigned trg_lookup = is_data_ ? run_ : year_;
    bool mu_trg=false;
    if (year_ == 2016) {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
      auto const& trg_objs_tk = event->GetPtrVec<TriggerObject>("triggerObjects_IsoTkMu24");
      mu_trg =
          IsFilterMatchedDR(m, trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3) ||
          IsFilterMatchedDR(m, trg_objs_tk, filters_IsoTkMu24_.Lookup(trg_lookup), 0.3);
    } else if (year_ == 2017){
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu27");
      mu_trg = IsFilterMatchedDR(m, trg_objs, filters_IsoMu27_.Lookup(trg_lookup), 0.3);
    } else {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
      mu_trg = IsFilterMatchedDR(m, trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3);
    }
    return mu_trg; 
  }

}
