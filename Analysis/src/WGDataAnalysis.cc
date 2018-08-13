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
#include "Acorn/NTupler/interface/Reduction.h"

namespace ac {

WGDataAnalysis::WGDataAnalysis(std::string const& name)
    : ModuleBase(name), fs_(nullptr), year_(2016), is_data_(true) {}

WGDataAnalysis::~WGDataAnalysis() { ; }

int WGDataAnalysis::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("WGDataAnalysis", "WGDataAnalysis");
    tree_->Branch("gen_proc", &gen_proc_);
    tree_->Branch("n_m", &n_m_);
    tree_->Branch("m0_pt", &m0_pt_);
    tree_->Branch("m0_eta", &m0_eta_);
    tree_->Branch("m0_phi", &m0_phi_);
    tree_->Branch("m0_trg", &m0_trg_);
    tree_->Branch("m1_pt", &m1_pt_);
    tree_->Branch("m1_eta", &m1_eta_);
    tree_->Branch("m1_phi", &m1_phi_);
    tree_->Branch("m0m1_M", &m0m1_M_);
    tree_->Branch("m0m1_dr", &m0m1_dr_);
    tree_->Branch("m0m1_os", &m0m1_os_);
    tree_->Branch("n_p", &n_p_);
    tree_->Branch("p0_pt", &p0_pt_);
    tree_->Branch("p0_eta", &p0_eta_);
    tree_->Branch("p0_phi", &p0_phi_);
    tree_->Branch("p0_chiso", &p0_chiso_);
    tree_->Branch("p0_medium", &p0_medium_);
    tree_->Branch("p0_tight", &p0_tight_);
    tree_->Branch("met", &met_);
    tree_->Branch("met_phi", &met_phi_);
    tree_->Branch("m0met_mt", &m0met_mt_);
    tree_->Branch("m0p0_dr", &m0p0_dr_);
    tree_->Branch("m0p0_dphi", &m0p0_dphi_);
    tree_->Branch("m0p0_M", &m0p0_M_);
    tree_->Branch("n_vm", &n_vm_);
    tree_->Branch("vm_p0_dr", &vm_p0_dr_);
    tree_->Branch("wt_def", &wt_def_);
    tree_->Branch("wt_pu", &wt_pu_);
    tree_->Branch("wt_m0", &wt_m0_);
    tree_->Branch("wt_trg_m0", &wt_trg_m0_);
    tree_->Branch("wt_m1", &wt_m1_);
    tree_->Branch("wt_p0", &wt_p0_);
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
  fns_["p_id_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_id_ratio")->functor(ws_->argSet("p_pt,p_eta")));
  return 0;
  }

  int WGDataAnalysis::Execute(TreeEvent* event) {

    SetDefaults();

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    auto muons = event->GetPtrVec<ac::Muon>("muons");
    auto photons = event->GetPtrVec<ac::Photon>("photons");
    auto mets = event->GetPtrVec<ac::Met>("pfType1Met");


    // Sub-process classification - maybe move this into a separate module
    if (gen_classify_ == "DY") {
      auto lhe_parts = event->GetPtrVec<ac::GenParticle>("lheParticles");
      unsigned n_ele = std::count_if(lhe_parts.begin(), lhe_parts.end(), [](ac::GenParticle *p) {
        return std::abs(p->pdgId()) == 11;
      });
      unsigned n_muo = std::count_if(lhe_parts.begin(), lhe_parts.end(), [](ac::GenParticle *p) {
        return std::abs(p->pdgId()) == 13;
      });
      unsigned n_tau = std::count_if(lhe_parts.begin(), lhe_parts.end(), [](ac::GenParticle *p) {
        return std::abs(p->pdgId()) == 15;
      });
      if (n_ele == 2) {
        gen_proc_ = 1;
      }
      if (n_muo == 2) {
        gen_proc_ = 2;
      }
      if (n_tau == 2) {
        gen_proc_ = 3;
      }
    }

    auto veto_muons = muons;
    ac::keep_if(veto_muons, [](ac::Muon const* m) {
      return m->pt() > 10. && fabs(m->eta()) < 2.4;
    });

    ac::keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && MuonPFIso(m) < 0.15;
    });

    // At this stage apply the medium Photon ID without the charged is cut
    ac::keep_if(photons, [](ac::Photon const* p) {
      return p->pt() > 30. && fabs(p->scEta()) < 2.5 &&
             (fabs(p->scEta()) < 1.4442 || fabs(p->scEta()) > 1.566) &&
             ac::PhotonIDIso(p, 2016, 1, false) && !p->hasPixelSeed();
    });

    boost::range::sort(muons, DescendingPt);
    boost::range::sort(photons, DescendingPt);

    if (muons.size() == 0) {
      return 1;
    }
    ac::Muon* m0 = muons[0];

    ac::keep_if(veto_muons, [&](ac::Muon *m) {
      return m != m0;
    });

    n_m_ = muons.size();
    n_vm_ = veto_muons.size();

    m0_pt_ = m0->pt();
    m0_pt_ = reduceMantissaToNbits(m0_pt_, 12);

    m0_eta_ = m0->eta();
    m0_phi_ = m0->phi();

    ac::Met* met = mets.at(0);
    met_ = met->pt();
    met_phi_ = met->phi();
    m0met_mt_ = ac::MT(m0, met);
    m0met_mt_ = reduceMantissaToNbits(m0met_mt_, 12);

    unsigned trg_lookup = is_data_ ? info->run() : year_;
    if (year_ == 2016) {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
      auto const& trg_objs_tk = event->GetPtrVec<TriggerObject>("triggerObjects_IsoTkMu24");
      m0_trg_ =
          IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3) ||
          IsFilterMatchedDR(muons[0], trg_objs_tk, filters_IsoTkMu24_.Lookup(trg_lookup), 0.3);
    } else {
      auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu27");
      m0_trg_ = IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu27_.Lookup(trg_lookup), 0.3);
    }

    if (muons.size() >= 2) {
      ac::Muon* m1 = muons[1];
      m1_pt_ = m1->pt();
      m1_eta_ = m1->eta();
      m1_phi_ = m1->phi();

      m0m1_M_ = (m0->vector() + m1->vector()).M();
      m0m1_dr_ = DeltaR(m0, m1);
      m0m1_os_ = m0->charge() != m1->charge();
    }

    n_p_ = photons.size();

    if (n_p_ >= 1) {
      ac::Photon* p0 = photons[0];
      p0_pt_ = p0->pt();
      p0_eta_ = p0->eta();
      p0_phi_ = p0->phi();
      p0_chiso_ = p0->chargedIso();
      p0_medium_ = p0->isMediumIdPhoton();
      p0_tight_ = p0->isTightIdPhoton();

      m0p0_dr_ =  ac::DeltaR(m0, p0);
      m0p0_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(m0->vector(), p0->vector());
      m0p0_M_ = (m0->vector() + p0->vector()).M();

      if (n_vm_ >= 1 && n_p_ >= 1) {
        boost::range::sort(veto_muons, [&](ac::Muon *x1, ac::Muon *x2) {
          return DeltaR(x1, p0) < DeltaR(x2, p0);
        });

        vm_p0_dr_ = ac::DeltaR(veto_muons[0], p0);
      }
    }

    if (!is_data_) {
      wt_def_ = info->nominalGenWeight() >= 0. ? +1. : -1.; // temporary until new ntuples fix EventInfo bug
      auto const& pu_info = event->GetPtrVec<PileupInfo>("pileupInfo");
      for (PileupInfo const* pu : pu_info) {
        if (pu->bunchCrossing() == 0) {
          wt_pu_ = RooFunc(fns_["pileup_ratio"], {pu->trueNumInteractions()});
          break;
        }
      }
      wt_m0_ = RooFunc(fns_["m_idisotrk_ratio"], {m0_pt_, m0_eta_});
      wt_trg_m0_ = RooFunc(fns_["m_trg_ratio"], {m0_pt_, m0_eta_});
      if (muons.size() >= 2) {
        wt_m1_ = RooFunc(fns_["m_idisotrk_ratio"], {m1_pt_, m1_eta_});
      }
      if (n_p_ >= 1) {
        wt_p0_ = RooFunc(fns_["p_id_ratio"], {p0_pt_, photons[0]->scEta()});
      }
    }

    // std::cout << "----\n";
    // auto genparts = event->GetPtrVec<ac::GenParticle>("genParticles");
    // for (auto p : genparts) {
    //   p->Print();
    // }
    tree_->Fill();



    return 0;
  }

  void WGDataAnalysis::SetDefaults() {
    gen_proc_ = 0;
    n_m_ = 0;
    m0_pt_ = 0.;
    m0_eta_ = 0.;
    m0_phi_ = 0.;
    m0_trg_ = false;
    m1_pt_ = 0.;
    m1_eta_ = 0.;
    m1_phi_ = 0.;
    m0m1_M_ = 0.;
    m0m1_dr_ = 0.;
    m0m1_os_ = false;
    n_p_ = 0;
    p0_pt_ = 0.;
    p0_eta_ = 0.;
    p0_phi_ = 0.;
    p0_chiso_ = 0.;
    p0_medium_ = false;
    p0_tight_ = false;
    met_ = 0.;
    met_phi_ = 0.;
    m0met_mt_ = 0.;
    m0p0_dr_ = 0.;
    m0p0_dphi_ = 0.;
    m0p0_M_ = 0.;
    n_vm_ = 0;
    vm_p0_dr_ = 0.;
    wt_def_ = 1.;
    wt_pu_ = 1.;
    wt_m0_ = 1.;
    wt_trg_m0_ = 1.;
    wt_m1_ = 1.;
    wt_p0_ = 1.;
  }

  int WGDataAnalysis::PostAnalysis() {
    return 0;
  }

  void WGDataAnalysis::PrintInfo() {}



}
