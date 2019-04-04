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
#include "Acorn/Analysis/interface/WGAnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/Reduction.h"

namespace ac {

WGDataAnalysis::WGDataAnalysis(std::string const& name)
    : ModuleBase(name),
      fs_(nullptr),
      year_(2016),
      is_data_(true),
      do_wg_gen_vars_(false),
      check_is_zg_(false),
      do_presel_(true) {}

WGDataAnalysis::~WGDataAnalysis() { ; }

int WGDataAnalysis::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("WGDataAnalysis", "WGDataAnalysis");
    tree_->Branch("run", &run_);
    tree_->Branch("gen_proc", &gen_proc_);
    tree_->Branch("gen_is_zg", &gen_is_zg_);
    tree_->Branch("n_vtx", &n_vtx_);
    tree_->Branch("n_pre_m", &n_pre_m_);
    tree_->Branch("n_pre_e", &n_pre_e_);
    tree_->Branch("n_veto_m", &n_veto_m_);
    tree_->Branch("n_veto_e", &n_veto_e_);
    tree_->Branch("l0_pt", &l0_pt_);
    tree_->Branch("l0_eta", &l0_eta_);
    tree_->Branch("l0_phi", &l0_phi_);
    tree_->Branch("l0_iso", &l0_iso_);
    tree_->Branch("l0_tkiso", &l0_tkiso_);
    tree_->Branch("l0_tight", &l0_tight_);
    tree_->Branch("l0_trg", &l0_trg_);
    tree_->Branch("l0_trg_2", &l0_trg_2_);
    tree_->Branch("l0_q", &l0_q_);
    tree_->Branch("l0_pdgid", &l0_pdgid_);
    tree_->Branch("l1_pt", &l1_pt_);
    tree_->Branch("l1_eta", &l1_eta_);
    tree_->Branch("l1_phi", &l1_phi_);
    tree_->Branch("l1_iso", &l1_iso_);
    tree_->Branch("l0l1_M", &l0l1_M_);
    tree_->Branch("l0l1_dr", &l0l1_dr_);
    tree_->Branch("l0l1_os", &l0l1_os_);
    tree_->Branch("n_pre_p", &n_pre_p_);
    tree_->Branch("p0_pt", &p0_pt_);
    tree_->Branch("p0_eta", &p0_eta_);
    tree_->Branch("p0_phi", &p0_phi_);
    tree_->Branch("p0_chiso", &p0_chiso_);
    tree_->Branch("p0_neiso", &p0_neiso_);
    tree_->Branch("p0_phiso", &p0_phiso_);
    tree_->Branch("p0_hovere", &p0_hovere_);
    tree_->Branch("p0_sigma", &p0_sigma_);
    tree_->Branch("p0_haspix", &p0_haspix_);
    tree_->Branch("p0_eveto", &p0_eveto_);
    tree_->Branch("p0_medium_noch", &p0_medium_noch_);
    tree_->Branch("p0_medium", &p0_medium_);
    tree_->Branch("p0_tight", &p0_tight_);
    tree_->Branch("p0_fsr", &p0_fsr_);
    tree_->Branch("p0_truth", &p0_truth_);
    tree_->Branch("met", &met_);
    tree_->Branch("met_phi", &met_phi_);
    tree_->Branch("xy_met", &xy_met_);
    tree_->Branch("xy_met_phi", &xy_met_phi_);
    tree_->Branch("puppi_met", &puppi_met_);
    tree_->Branch("puppi_met_phi", &puppi_met_phi_);
    tree_->Branch("l0met_mt", &l0met_mt_);
    tree_->Branch("l0p0_dr", &l0p0_dr_);
    tree_->Branch("l0p0_dphi", &l0p0_dphi_);
    tree_->Branch("l0p0_M", &l0p0_M_);
    tree_->Branch("reco_phi", &reco_phi_);
    tree_->Branch("reco_sphi", &reco_sphi_);
    tree_->Branch("reco_xy_phi", &reco_xy_phi_);
    tree_->Branch("reco_xy_sphi", &reco_xy_sphi_);
    tree_->Branch("reco_puppi_phi", &reco_puppi_phi_);
    tree_->Branch("reco_puppi_sphi", &reco_puppi_sphi_);
    tree_->Branch("n_vm", &n_vm_);
    tree_->Branch("vm_p0_dr", &vm_p0_dr_);
    tree_->Branch("wt_def", &wt_def_);
    tree_->Branch("wt_pf", &wt_pf_);
    tree_->Branch("wt_pu", &wt_pu_);
    tree_->Branch("wt_l0", &wt_l0_);
    tree_->Branch("wt_trg_l0", &wt_trg_l0_);
    tree_->Branch("wt_l1", &wt_l1_);
    tree_->Branch("wt_p0", &wt_p0_);
    tree_->Branch("wt_p0_fake", &wt_p0_fake_);
    tree_->Branch("wt_p0_e_fake", &wt_p0_e_fake_);
    tree_->Branch("is_wg_gen", &is_wg_gen_);
    tree_->Branch("gen_nparts", &gen_nparts_);
    tree_->Branch("gen_pdgid", &gen_pdgid_);
    tree_->Branch("gen_l0_match", &gen_l0_match_);
    tree_->Branch("gen_p0_match", &gen_p0_match_);
    tree_->Branch("lhe_l0_pt", &lhe_l0_pt_);
    tree_->Branch("lhe_l0_eta", &lhe_l0_eta_);
    tree_->Branch("lhe_p0_pt", &lhe_p0_pt_);
    tree_->Branch("lhe_p0_eta", &lhe_p0_eta_);
    tree_->Branch("gen_p0_pt", &gen_p0_pt_);
    tree_->Branch("gen_p0_eta", &gen_p0_eta_);
    tree_->Branch("gen_phi", &gen_phi_);
    tree_->Branch("gen_sphi", &gen_sphi_);
    tree_->Branch("true_phi", &true_phi_);
    tree_->Branch("gen_l0_q", &gen_l0_q_);
    tree_->Branch("gen_l0_pt", &gen_l0_pt_);
    tree_->Branch("gen_l0_eta", &gen_l0_eta_);
    tree_->Branch("gen_met", &gen_met_);
    tree_->Branch("gen_l0p0_dr", &gen_l0p0_dr_);
    // tree_->Branch("gen_dy_mll", &gen_dy_mll_);
    // tree_->Branch("gen_n2_pt", &gen_n2_pt_);
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

    filters_Mu50_ =
        LookupFilter({{272023, "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"}});

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

    filters_Mu50_ =
        LookupFilter({{2016, "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"}});

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
  fns_["m_idisotrk_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_idisotrk_ratio")->functor(ws_->argSet("m_pt,m_eta")));
  fns_["m_trg_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_trg_ratio")->functor(ws_->argSet("m_pt,m_eta")));
  fns_["e_gsfidiso_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_gsfidiso_ratio")->functor(ws_->argSet("e_pt,e_eta")));
  fns_["e_trg_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_trg_ratio")->functor(ws_->argSet("e_pt,e_eta")));
  fns_["p_id_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_id_ratio")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["e_p_fake_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_p_fake_ratio")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_fake_ratio_m_chn"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_fake_ratio_m_chn")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_fake_ratio_e_chn"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_fake_ratio_e_chn")->functor(ws_->argSet("p_pt,p_eta")));
  return 0;
  }

  int WGDataAnalysis::Execute(TreeEvent* event) {

    SetDefaults();

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    run_ = info->run();
    n_vtx_ = info->numVertices();
    metfilters_ = info->metfilters().to_ulong();

    auto muons = event->GetPtrVec<ac::Muon>("muons");
    auto electrons = event->GetPtrVec<ac::Electron>("electrons");
    auto photons = event->GetPtrVec<ac::Photon>("photons");

    // TODO: temporary hack
    for (Photon *p : photons) {
      PhotonIsoCorrector(p, n_vtx_);
    }
    auto mets = event->GetPtrVec<ac::Met>("pfType1Met");
    auto puppi_mets = event->GetPtrVec<ac::Met>("puppiMet");

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

    auto pre_muons = ac::copy_keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && fabs(m->dxy()) < 0.05 &&
             fabs(m->dz()) < 0.2;
    });

    // TODO: update with SC eta
    auto pre_electrons = ac::copy_keep_if(electrons, [](ac::Electron const* e) {
      return e->pt() > 35. && fabs(e->eta()) < 2.5 && e->isCutBasedMediumElectron() &&
             (fabs(e->eta()) < 1.4442 || fabs(e->eta()) > 1.566) && ElectronIPCuts(e);
    });

    // At this stage apply the medium Photon ID without the charged iso cut
    auto pre_photons = ac::copy_keep_if(photons, [&](ac::Photon const* p) {
      return p->pt() > 30. && fabs(p->scEta()) < 2.5 &&
             (fabs(p->scEta()) < 1.4442 || fabs(p->scEta()) > 1.566) && PhotonIDIso(p, year_, 1, false, false);
    });

    // TODO: update with proper loose ID
    auto veto_muons = ac::copy_keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 10. && fabs(m->eta()) < 2.4 && fabs(m->dxy()) < 0.05 &&
             fabs(m->dz()) < 0.2;
    });

    auto veto_electrons = ac::copy_keep_if(electrons, [](ac::Electron const* e) {
      return e->pt() > 10. && fabs(e->eta()) < 2.5 && e->isCutBasedLooseElectron() &&
             (fabs(e->eta()) < 1.4442 || fabs(e->eta()) > 1.566) && ElectronIPCuts(e);
    });

    boost::range::sort(pre_muons, DescendingPt);
    boost::range::sort(pre_electrons, DescendingPt);
    boost::range::sort(pre_photons, DescendingPt);

    // Resolve the case where we have at n_m >= 1 and n_e >= 1
    ac::Muon* m0 = pre_muons.size() ? pre_muons[0] : nullptr;
    ac::Electron* e0 = pre_electrons.size() ? pre_electrons[0] : nullptr;

    ac::Candidate *l0 = nullptr;
    if (m0 && e0) {
      if (m0->pt() > e0->pt()) {
        l0 = m0;
        e0 = nullptr;
      } else {
        l0 = e0;
        m0 = nullptr;
      }
    } else if (m0) {
      l0 = m0;
    } else if (e0) {
      l0 = e0;
    }

    if (m0) {
      l0_pdgid_ = 13;
    } else if (e0) {
      l0_pdgid_ = 11;
    }

    if (l0) {
      ac::keep_if(pre_photons, [&](ac::Photon const* p) {
        return DeltaR(p, l0) > 0.3;
      });
    }

    ac::Photon* p0 = pre_photons.size() ? pre_photons[0] : nullptr;

    // Here are the conditions for skipping this event and not writing anything to the tree
    if (do_presel_ && (l0 == nullptr || p0 == nullptr)) {
      return 1;
    }

    n_pre_m_ = pre_muons.size();
    n_pre_e_ = pre_electrons.size();
    n_pre_p_ = pre_photons.size();

    n_veto_m_ = veto_muons.size();
    n_veto_e_ = veto_electrons.size();

    // Always remove m0 from the list of veto muons
    if (m0) {
      ac::keep_if(veto_muons, [&](ac::Muon *m) {
        return m != m0;
      });
    }

    n_vm_ = veto_muons.size();

    // Calculate the truth vars for W+g events
    if (do_wg_gen_vars_) {
      auto gen_parts = event->GetPtrVec<ac::GenParticle>("genParticles");
      auto lhe_parts = event->GetPtrVec<ac::GenParticle>("lheParticles");
      auto gen_met = event->GetPtrVec<ac::Met>("genMet")[0];
      WGGenParticles parts = ProduceWGGenParticles(lhe_parts, gen_parts);
      gen_nparts_ = parts.nparts;
      if (parts.lhe_lep) {
        gen_pdgid_ = std::abs(parts.lhe_lep->pdgId());
        lhe_l0_pt_ = parts.lhe_lep->pt();
        lhe_l0_eta_ = parts.lhe_lep->eta();
      }
      if (parts.lhe_pho) {
        lhe_p0_pt_ = parts.lhe_pho->pt();
        lhe_p0_eta_ = parts.lhe_pho->eta();
      }
      if (parts.ok) {
        is_wg_gen_ = true;
        WGSystem gen_sys = ProduceWGSystem(*parts.gen_lep, *gen_met, *parts.gen_pho, true, rng, false);
        WGSystem gen_true_sys = ProduceWGSystem(*parts.gen_lep, *parts.gen_neu, *parts.gen_pho, false, rng, false);
        gen_p0_pt_ = parts.gen_pho->pt();
        gen_p0_eta_ = parts.gen_pho->eta();
        gen_l0_q_ = parts.gen_lep->charge();
        gen_phi_ = reduceMantissaToNbits(gen_sys.Phi(parts.gen_lep->charge() > 0), 12);
        gen_sphi_ = reduceMantissaToNbits(gen_sys.SymPhi(parts.gen_lep->charge() > 0), 12);
        gen_l0_pt_ = parts.gen_lep->pt();
        gen_l0_eta_ = parts.gen_lep->eta();
        gen_met_ = gen_met->pt();
        gen_l0p0_dr_ = reduceMantissaToNbits(ac::DeltaR(parts.gen_lep, parts.gen_pho), 12);
        true_phi_ = reduceMantissaToNbits(gen_true_sys.Phi(parts.gen_lep->charge() > 0), 12);

        // Now try and match to the reco objects, if they exist
        if (l0 && ac::DeltaR(l0, parts.gen_lep) < 0.3) {
          gen_l0_match_ = true;
        }
        if (p0 && ac::DeltaR(p0, parts.gen_pho) < 0.3) {
          gen_p0_match_ = true;
        }
      }
    }

    ac::Met* met = mets.at(0);
    met_ = met->pt();
    met_phi_ = met->phi();

    ac::Met* xy_met = mets.at(1);
    xy_met_ = xy_met->pt();
    xy_met_phi_ = xy_met->phi();

    ac::Met* puppi_met = puppi_mets.at(0);
    puppi_met_ = puppi_met->pt();
    puppi_met_phi_ = puppi_met->phi();

    if (l0) {
      l0_pt_ = l0->pt();
      l0_eta_ = l0->eta();
      l0_phi_ = l0->phi();
      l0_q_ = l0->charge();
      l0met_mt_ = ac::MT(l0, met);
      l0met_mt_ = reduceMantissaToNbits(l0met_mt_, 12);

      if (m0) {
        l0_iso_ = MuonPFIso(m0);
        l0_tkiso_ = m0->pfIsoSumChargedHadronPt() / m0->pt();
        l0_tight_ = m0->isTightMuon();
      } else if (e0) {
        l0_iso_ = e0->relativeEAIso();
        l0_tight_ = e0->isCutBasedTightElectron();
      }

      unsigned trg_lookup = is_data_ ? info->run() : year_;
      // Muon trigger strategies:
      // 1) IsoMu24/TkMu24/IsoMu27
      // 2) (1) || Mu50/TkMu50 if pT >~ 55 GeV
      if (m0) {
        if (year_ == 2016) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
          auto const& trg_objs_tk = event->GetPtrVec<TriggerObject>("triggerObjects_IsoTkMu24");
          l0_trg_ =
              IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3) ||
              IsFilterMatchedDR(muons[0], trg_objs_tk, filters_IsoTkMu24_.Lookup(trg_lookup), 0.3);
        } else if (year_ == 2017) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu27");
          l0_trg_ = IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu27_.Lookup(trg_lookup), 0.3);
        } else if (year_ == 2018) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
          l0_trg_ = IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3);
        }
        auto const& trg_objs_Mu50 = event->GetPtrVec<TriggerObject>("triggerObjects_Mu50");
        l0_trg_2_ = IsFilterMatchedDR(muons[0], trg_objs_Mu50, filters_Mu50_.Lookup(trg_lookup), 0.3);
      }
      if (e0) {
        if (year_ == 2016) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele27_WPTight_Gsf");
          l0_trg_ = IsFilterMatchedDR(e0, trg_objs, filters_Ele27_.Lookup(trg_lookup), 0.3);
        } else if (year_ == 2017) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf");
          l0_trg_ = IsFilterMatchedDR(e0, trg_objs, filters_Ele32_.Lookup(trg_lookup), 0.3);
          if (is_data_ && info->run() < 302026) {
            auto const& trg_objs_alt = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf_L1DoubleEG");
            bool l0_trg_alt = IsFilterMatchedDR(e0, trg_objs_alt, filters_Ele32_L1DoubleEG_.Lookup(trg_lookup), 0.3) &&
                              IsFilterMatchedDR(e0, trg_objs_alt, filters_Ele32_L1DoubleEG_seed_.Lookup(trg_lookup), 0.3);
            l0_trg_ = l0_trg_ || l0_trg_alt;
          }
        } else if (year_ == 2018) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf");
          l0_trg_ = IsFilterMatchedDR(e0, trg_objs, filters_Ele32_.Lookup(trg_lookup), 0.3);
        }
      }
    }

    if (m0 && muons.size() >= 2) {
      ac::Muon* m1 = muons[1];
      l1_pt_ = m1->pt();
      l1_eta_ = m1->eta();
      l1_phi_ = m1->phi();
      l1_iso_ = MuonPFIso(m1);

      l0l1_M_ = (m0->vector() + m1->vector()).M();
      l0l1_dr_ = DeltaR(m0, m1);
      l0l1_os_ = m0->charge() != m1->charge();
    } else if (e0 && electrons.size() >= 2) {
      ac::Electron* e1 = electrons[1];
      l1_pt_ = e1->pt();
      l1_eta_ = e1->eta();
      l1_phi_ = e1->phi();
      l1_iso_ = e1->relativeEAIso();

      l0l1_M_ = (e0->vector() + e1->vector()).M();
      l0l1_dr_ = DeltaR(e0, e1);
      l0l1_os_ = e0->charge() != e1->charge();
    }
    l0l1_M_ = reduceMantissaToNbits(l0l1_M_, 12);
    l0l1_dr_ = reduceMantissaToNbits(l0l1_dr_, 12);

    if (p0) {
      p0_pt_ = p0->pt();
      p0_eta_ = p0->eta();
      p0_phi_ = p0->phi();
      p0_chiso_ = p0->chargedIso();
      p0_neiso_ = p0->neutralHadronIso();
      p0_phiso_ = p0->photonIso();
      p0_hovere_ = p0->hadTowOverEm();
      p0_sigma_ = p0->full5x5SigmaIetaIeta();
      p0_haspix_ = p0->hasPixelSeed();
      p0_eveto_ = p0->passElectronVeto();
      p0_medium_noch_ = ac::PhotonIDIso(p0, year_, 1, false, false);
      p0_medium_ = p0->isMediumIdPhoton();
      p0_tight_ = p0->isTightIdPhoton();

      if (n_vm_ >= 1) {
        boost::range::sort(veto_muons, [&](ac::Muon* x1, ac::Muon* x2) {
          return DeltaR(x1, p0) < DeltaR(x2, p0);
        });
        vm_p0_dr_ = ac::DeltaR(veto_muons[0], p0);
      }
    }

    if (p0 && l0) {
      l0p0_dr_ =  reduceMantissaToNbits(ac::DeltaR(l0, p0), 12);
      l0p0_dphi_ = reduceMantissaToNbits(ROOT::Math::VectorUtil::DeltaPhi(l0->vector(), p0->vector()), 12);
      l0p0_M_ = reduceMantissaToNbits((l0->vector() + p0->vector()).M(), 12);

      WGSystem reco_sys = ProduceWGSystem(*l0, *met, *p0, true, rng, false);
      reco_phi_ = reduceMantissaToNbits(reco_sys.Phi(l0->charge()), 7);
      reco_sphi_ = reduceMantissaToNbits(reco_sys.SymPhi(l0->charge()), 7);

      WGSystem reco_xy_sys = ProduceWGSystem(*l0, *xy_met, *p0, true, rng, false);
      reco_xy_phi_ = reduceMantissaToNbits(reco_xy_sys.Phi(l0->charge()), 7);
      reco_xy_sphi_ = reduceMantissaToNbits(reco_xy_sys.SymPhi(l0->charge()), 7);

      WGSystem reco_puppi_sys = ProduceWGSystem(*l0, *puppi_met, *p0, true, rng, false);
      reco_puppi_phi_ = reduceMantissaToNbits(reco_puppi_sys.Phi(l0->charge()), 7);
      reco_puppi_sphi_ = reduceMantissaToNbits(reco_puppi_sys.SymPhi(l0->charge()), 7);

      if (m0) {
        wt_p0_fake_ = RooFunc(fns_["p_fake_ratio_m_chn"], {p0_pt_, p0->scEta()});
      }
      if (e0) {
        wt_p0_fake_ = RooFunc(fns_["p_fake_ratio_e_chn"], {p0_pt_, p0->scEta()});
      }
    }

    if (!is_data_) {
      // wt_def_ = info->nominalGenWeight() >= 0. ? +1. : -1.; // temporary until new ntuples fix EventInfo bug
      wt_def_ = info->totalWeight();
      auto const& pu_info = event->GetPtrVec<PileupInfo>("pileupInfo");
      for (PileupInfo const* pu : pu_info) {
        if (pu->bunchCrossing() == 0) {
          wt_pu_ = RooFunc(fns_["pileup_ratio"], {pu->trueNumInteractions()});
          break;
        }
      }
      if (info->userDoubles().size() >= 1) {
        wt_pf_ = info->userDoubles().at(0);
      }
      if (m0) {
        wt_l0_ = RooFunc(fns_["m_idisotrk_ratio"], {l0_pt_, l0_eta_});
        wt_trg_l0_ = RooFunc(fns_["m_trg_ratio"], {l0_pt_, l0_eta_});
      }
      if (e0) {
        wt_l0_ = RooFunc(fns_["e_gsfidiso_ratio"], {l0_pt_, l0_eta_});
        wt_trg_l0_ = RooFunc(fns_["e_trg_ratio"], {l0_pt_, l0_eta_});
      }
      if (m0 && muons.size() >= 2) {
        wt_l1_ = RooFunc(fns_["m_idisotrk_ratio"], {l1_pt_, l1_eta_});
      }
      auto genparts = event->GetPtrVec<ac::GenParticle>("genParticles");
      if (p0) {
        auto gen_st1_dr = ac::copy_keep_if(genparts, [&](ac::GenParticle *p) {
          return p->status() == 1 && DeltaR(p, p0) < 0.3;
        });
        auto gen_photons = ac::copy_keep_if(gen_st1_dr, [&](ac::GenParticle *p) {
          return p->pt() > 10. && p->pdgId() == 22;
        });
        auto prompt_gen_photons = ac::copy_keep_if(gen_photons, [&](ac::GenParticle *p) {
          return p->statusFlags().isPrompt();
        });
        auto hadronic_gen_photons = ac::copy_keep_if(gen_photons, [&](ac::GenParticle *p) {
          return p->statusFlags().isDirectHadronDecayProduct();
        });
        auto all_prompt_gen_elecs = ac::copy_keep_if(gen_st1_dr, [&](ac::GenParticle *p) {
          return p->pt() > 0. && std::abs(p->pdgId()) == 11 && p->statusFlags().isPrompt();
        });
        auto all_prompt_gen_muons = ac::copy_keep_if(gen_st1_dr, [&](ac::GenParticle *p) {
          return p->pt() > 0. && std::abs(p->pdgId()) == 13 && p->statusFlags().isPrompt();
        });
        auto prompt_gen_elecs = ac::copy_keep_if(all_prompt_gen_elecs, [&](ac::GenParticle *p) {
          return p->pt() > 10.;
        });
        auto prompt_gen_muons = ac::copy_keep_if(all_prompt_gen_muons, [&](ac::GenParticle *p) {
          return p->pt() > 10.;
        });

        // Only prompt e/mu, pT > 10 -> lep->photon fake
        if (prompt_gen_elecs.size() > 0 && prompt_gen_photons.size() == 0) {
          p0_truth_ = 2;
        } else if (prompt_gen_muons.size() > 0 && prompt_gen_photons.size() == 0) {
          p0_truth_ = 3;
        // Prompt photon + any lepton -> FSR photon
        } else if (all_prompt_gen_elecs.size() > 0 && prompt_gen_photons.size() > 0) {
          p0_truth_ = 4;
        } else if (all_prompt_gen_muons.size() > 0 && prompt_gen_photons.size() > 0) {
          p0_truth_ = 5;
        } else if (hadronic_gen_photons.size() > 0) {
          p0_truth_ = 6;
        } else if (prompt_gen_photons.size() > 0) {
          p0_truth_ = 1;
        }
        for (auto const& p_pho : prompt_gen_photons) {
          for (auto const& p : genparts) {
            if (ac::IsChargedLepton(*p) && ac::contains(p->daughters(), p_pho->index())) {
              p0_fsr_ = true;
              break;
            }
          }
          if (p0_fsr_) break;
        }
        if (ac::contains({1, 2, 4, 5}, p0_truth_)) {
          wt_p0_ = RooFunc(fns_["p_id_ratio"], {p0_pt_, p0->scEta()});
        }
        if (ac::contains({2}, p0_truth_)) {
          wt_p0_e_fake_ = RooFunc(fns_["e_p_fake_ratio"], {p0_pt_, p0->eta()});
        }
      }

      // Check if this is the ZG sample phasespace
      if (check_is_zg_) {
        ROOT::Math::PtEtaPhiMVector dy_parts;
        std::vector<double> pts;
        for (auto const& p : genparts) {
          if (p->statusFlags().isPrompt() && IsChargedLepton(*p) && p->statusFlags().isLastCopy()) {
            dy_parts += p->vector();
            pts.push_back(p->pt());
          }
        }
        std::sort(pts.begin(), pts.end(), [](double d1, double d2) {return d1 > d2; });
        if (pts.size() >= 2 && pts.at(1) > 15. && dy_parts.M() > 30.) {
          gen_is_zg_ = true;
        }
      }
    }
    tree_->Fill();
    return 0;
  }

  void WGDataAnalysis::SetDefaults() {
    run_ = 0;
    gen_proc_ = 0;
    gen_is_zg_ = false;
    n_vtx_ = 0;
    metfilters_ = 0;
    n_pre_m_ = 0;
    n_pre_e_ = 0;
    n_veto_m_ = 0;
    n_veto_e_ = 0;
    l0_pt_ = 0.;
    l0_eta_ = 0.;
    l0_phi_ = 0.;
    l0_iso_ = 0.;
    l0_tkiso_ = 0.;
    l0_tight_ = false;
    l0_trg_ = false;
    l0_trg_2_ = false;
    l0_pdgid_ = 0;
    l1_pt_ = 0.;
    l1_eta_ = 0.;
    l1_phi_ = 0.;
    l1_iso_ = 0.;
    l0l1_M_ = 0.;
    l0l1_dr_ = 0.;
    l0l1_os_ = false;
    n_pre_p_ = 0;
    p0_pt_ = 0.;
    p0_eta_ = 0.;
    p0_phi_ = 0.;
    p0_chiso_ = 0.;
    p0_neiso_ = 0.;
    p0_phiso_ = 0.;
    p0_hovere_ = 0.;
    p0_sigma_ = 0.;
    p0_haspix_ = false;
    p0_eveto_ = false;
    p0_medium_noch_ = false;
    p0_medium_ = false;
    p0_fsr_ = false;
    p0_tight_ = false;
    p0_truth_ = 0;
    met_ = 0.;
    met_phi_ = 0.;
    xy_met_ = 0.;
    xy_met_phi_ = 0.;
    puppi_met_ = 0.;
    puppi_met_phi_ = 0.;
    l0met_mt_ = 0.;
    l0p0_dr_ = 0.;
    l0p0_dphi_ = 0.;
    l0p0_M_ = 0.;
    reco_phi_ = 0.;
    reco_sphi_ = 0.;
    reco_xy_phi_ = 0.;
    reco_xy_sphi_ = 0.;
    reco_puppi_phi_ = 0.;
    reco_puppi_sphi_ = 0.;
    n_vm_ = 0;
    vm_p0_dr_ = 0.;
    wt_def_ = 1.;
    wt_pf_ = 1.;
    wt_pu_ = 1.;
    wt_l0_ = 1.;
    wt_trg_l0_ = 1.;
    wt_l1_ = 1.;
    wt_p0_ = 1.;
    wt_p0_fake_ = 1.;
    wt_p0_e_fake_ = 1.;
    is_wg_gen_ = false;
    gen_pdgid_ = 0;
    gen_nparts_ = 0;
    gen_l0_match_ = false;
    gen_p0_match_ = false;
    lhe_l0_pt_ = 0.;
    lhe_l0_eta_ = 0.;
    lhe_p0_pt_ = 0.;
    lhe_p0_eta_ = 0.;
    gen_p0_pt_ = 0.;
    gen_p0_eta_ = 0.;
    gen_phi_ = 0.;
    gen_sphi_ = 0.;
    gen_l0_q_ = 0;
    gen_l0_pt_ = 0.;
    gen_l0_eta_ = 0.;
    gen_met_ = 0.;
    gen_l0p0_dr_ = 0.;
    true_phi_ = 0.;
    // gen_dy_mll_ = 0.;
    // gen_n2_pt_ = 0.;
  }

  int WGDataAnalysis::PostAnalysis() {
    return 0;
  }

  void WGDataAnalysis::PrintInfo() {}

  void WGDataAnalysis::PhotonIsoCorrector(ac::Photon *p, unsigned nvertices) {
    double nv = double(nvertices);
    double rho = 0.291781 + nv * 0.736966 + nv * nv * -0.00152866;

    double a_eta = std::abs(p->scEta());
    double chiso = p->chargedIso();
    double neiso = p->neutralHadronIso();
    double phiso = p->photonIso();

    double ch_ea = 0.;
    double ne_ea = 0.;
    double ph_ea = 0.;

    if (a_eta < 1.0) {
      ch_ea = 0.0112;
      ne_ea = 0.0668;
      ph_ea = 0.1113;
    } else if (a_eta < 1.479) {
      ch_ea = 0.0108;
      ne_ea = 0.1054;
      ph_ea = 0.0953;
    } else if (a_eta < 2.0) {
      ch_ea = 0.0106;
      ne_ea = 0.0786;
      ph_ea = 0.0619;
    } else if (a_eta < 2.2) {
      ch_ea = 0.01002;
      ne_ea = 0.0233;
      ph_ea = 0.0837;
    } else if (a_eta < 2.3) {
      ch_ea = 0.0098;
      ne_ea = 0.0078;
      ph_ea = 0.1070;
    } else if (a_eta < 2.4) {
      ch_ea = 0.0089;
      ne_ea = 0.0028;
      ph_ea = 0.1212;
    } else {
      ch_ea = 0.0087;
      ne_ea = 0.0137;
      ph_ea = 0.1466;
    }
    chiso = std::max(chiso - ch_ea * rho, 0.);
    neiso = std::max(neiso - ne_ea * rho, 0.);
    phiso = std::max(phiso - ph_ea * rho, 0.);
    p->setChargedIso(chiso);
    p->setNeutralHadronIso(neiso);
    p->setPhotonIso(phiso);
  }
}
