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
#include "Acorn/NTupler/interface/PFJet.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/Reduction.h"
#include "Acorn/NTupler/src/CMS_2020_PAS_SMP_20_005.h"

namespace ac {

WGDataAnalysis::WGDataAnalysis(std::string const& name)
    : ModuleBase(name),
      fs_(nullptr),
      year_(2016),
      is_data_(true),
      do_wg_gen_vars_(false),
      check_is_zg_(false),
      check_is_wwg_(false),
      check_gen_mll_(false),
      do_presel_(true),
      only_wg_(false),
      var_set_(1),
      correct_e_energy_(-1),
      correct_p_energy_(-1),
      correct_m_energy_(-1),
      shift_met_(-1),
      scale_weights_(0),
      ps_weights_(0) {}

WGDataAnalysis::~WGDataAnalysis() { ; }

int WGDataAnalysis::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("WGDataAnalysis", "WGDataAnalysis");
    // The default AutoFlush (30MB) seems to cause large
    // memory usage at times (probably due to some branches
    // being basically empty). Reduce it by a factor of 10 here:
    tree_->SetAutoFlush(-6000000);
    // var sets:
    // 0 = essential
    // 1 = nominal
    // 2 = all

    if (var_set_ >= 0) {
      tree_->Branch("metfilters", &metfilters_);
      if (check_is_zg_) {
        tree_->Branch("gen_is_zg", &gen_is_zg_);
      }
      if (check_is_wwg_) {
        tree_->Branch("gen_proc", &gen_proc_);
      }
      tree_->Branch("n_vtx", &n_vtx_);

      tree_->Branch("n_pre_m", &n_pre_m_);
      tree_->Branch("n_pre_e", &n_pre_e_);
      tree_->Branch("n_veto_m", &n_veto_m_);
      tree_->Branch("n_veto_e", &n_veto_e_);
      tree_->Branch("n_alt_veto_m", &n_alt_veto_m_);
      tree_->Branch("n_alt_veto_e", &n_alt_veto_e_);
      tree_->Branch("n_qcd_j", &n_qcd_j_);

      tree_->Branch("l0_pt", &l0_pt_);
      tree_->Branch("l0_eta", &l0_eta_);
      tree_->Branch("l0_phi", &l0_phi_);
      tree_->Branch("l0_iso", &l0_iso_);
      tree_->Branch("l0_nominal", &l0_nominal_);
      tree_->Branch("l0_medium", &l0_medium_);
      tree_->Branch("l0_tight", &l0_tight_);
      tree_->Branch("l0_trg", &l0_trg_);
      tree_->Branch("l0_q", &l0_q_);
      tree_->Branch("l0_pdgid", &l0_pdgid_);

      tree_->Branch("l1_pt", &l1_pt_);
      // tree_->Branch("l1_iso", &l1_iso_);
      tree_->Branch("l0l1_M", &l0l1_M_);
      tree_->Branch("l0l1_dr", &l0l1_dr_);
      tree_->Branch("l0l1_pt", &l0l1_pt_);
      tree_->Branch("l0l1_os", &l0l1_os_);

      tree_->Branch("n_pre_p", &n_pre_p_);
      tree_->Branch("p0_pt", &p0_pt_);
      tree_->Branch("p0_eta", &p0_eta_);
      tree_->Branch("p0_phi", &p0_phi_);
      tree_->Branch("p0_chiso", &p0_chiso_);
      tree_->Branch("p0_sigma", &p0_sigma_);
      tree_->Branch("p0_haspix", &p0_haspix_);
      tree_->Branch("p0_eveto", &p0_eveto_);
      tree_->Branch("p0_medium_noch", &p0_medium_noch_);
      tree_->Branch("p0_loose", &p0_loose_);
      tree_->Branch("p0_medium", &p0_medium_);
      tree_->Branch("p0_truth", &p0_truth_);

      // tree_->Branch("j0_pt", &j0_pt_);
      // tree_->Branch("j0_eta", &j0_eta_);
      tree_->Branch("l0j0_dphi", &l0j0_dphi_);

      tree_->Branch("met", &met_);
      // tree_->Branch("tk_met", &tk_met_);
      // tree_->Branch("tk_met_phi", &tk_met_phi_);
      tree_->Branch("puppi_met", &puppi_met_);

      tree_->Branch("l0met_mt", &l0met_mt_);

      tree_->Branch("l0p0_dr", &l0p0_dr_);
      // tree_->Branch("l0p0_dphi", &l0p0_dphi_);
      tree_->Branch("l0p0_M", &l0p0_M_);
      tree_->Branch("mt_cluster", &mt_cluster_);

      // tree_->Branch("reco_phi", &reco_phi_);
      // tree_->Branch("reco_phi_f", &reco_phi_f_);
      // tree_->Branch("reco_tk_phi", &reco_tk_phi_);
      // tree_->Branch("reco_tk_phi_f", &reco_tk_phi_f_);
      // tree_->Branch("reco_puppi_phi", &reco_puppi_phi_);
      tree_->Branch("reco_puppi_phi_f", &reco_puppi_phi_f_);

      tree_->Branch("wt_def", &wt_def_);
      tree_->Branch("wt_pf", &wt_pf_);
      tree_->Branch("wt_pu", &wt_pu_);
      tree_->Branch("wt_l0", &wt_l0_);
      tree_->Branch("wt_trg_l0", &wt_trg_l0_);
      tree_->Branch("wt_l1", &wt_l1_);
      tree_->Branch("wt_p0", &wt_p0_);
      tree_->Branch("wt_p0_fake", &wt_p0_fake_);
      tree_->Branch("wt_p0_highpt_fake", &wt_p0_highpt_fake_);
      tree_->Branch("wt_p0_e_fake", &wt_p0_e_fake_);
      tree_->Branch("wt_l0_fake", &wt_l0_fake_);

      if (do_wg_gen_vars_) {
        tree_->Branch("gen_p0_pt", &gen_p0_pt_);
        tree_->Branch("gen_p0_eta", &gen_p0_eta_);
        tree_->Branch("gen_phi", &gen_phi_);
        tree_->Branch("gen_phi_f", &gen_phi_f_);
        tree_->Branch("true_phi", &true_phi_);
        tree_->Branch("true_phi_f", &true_phi_f_);

        tree_->Branch("gen_l0_q", &gen_l0_q_);
        tree_->Branch("gen_l0_pt", &gen_l0_pt_);
        tree_->Branch("gen_l0_eta", &gen_l0_eta_);
        tree_->Branch("gen_met", &gen_met_);
        tree_->Branch("gen_l0p0_dr", &gen_l0p0_dr_);
        tree_->Branch("gen_wg_M", &gen_wg_M_);
        tree_->Branch("gen_mt_cluster", &gen_mt_cluster_);

        tree_->Branch("lhe_frixione", &lhe_frixione_);
      }
    }

    if (var_set_ >= 1) {
      tree_->Branch("gen_mll", &gen_mll_);


      tree_->Branch("l1_eta", &l1_eta_);
      tree_->Branch("l1_phi", &l1_phi_);

      tree_->Branch("met_phi", &met_phi_);
      tree_->Branch("puppi_met_phi", &puppi_met_phi_);

      tree_->Branch("wt_sc_0", &wt_sc_0_);
      tree_->Branch("wt_sc_1", &wt_sc_1_);
      tree_->Branch("wt_sc_2", &wt_sc_2_);
      tree_->Branch("wt_sc_3", &wt_sc_3_);
      tree_->Branch("wt_sc_4", &wt_sc_4_);
      tree_->Branch("wt_sc_5", &wt_sc_5_);

      int npdf_tot = 0;
      for (unsigned ipdfset = 0; ipdfset < pdf_begin_.size(); ++ipdfset) {
        npdf_tot += (pdf_end_.at(ipdfset) - pdf_begin_.at(ipdfset) + 1);
      }
      wt_pdf_.resize(npdf_tot);

      int ipdf = 0;
      for (unsigned ipdfset = 0; ipdfset < pdf_begin_.size(); ++ipdfset) {
        int npdf = pdf_end_.at(ipdfset) - pdf_begin_.at(ipdfset) + 1;
        for (int ip = 0; ip < npdf; ++ip) {
          tree_->Branch(TString::Format("wt_pdf_%i", ipdf), &(wt_pdf_[ipdf]));
          ++ipdf;
        }
      }

      tree_->Branch("wt_isr_hi", &wt_isr_hi_);
      tree_->Branch("wt_isr_lo", &wt_isr_lo_);
      tree_->Branch("wt_fsr_hi", &wt_fsr_hi_);
      tree_->Branch("wt_fsr_lo", &wt_fsr_lo_);

      tree_->Branch("wt_pu_hi", &wt_pu_hi_);
      tree_->Branch("wt_pu_lo", &wt_pu_lo_);
      tree_->Branch("wt_pf_hi", &wt_pf_hi_);
      tree_->Branch("wt_pf_lo", &wt_pf_lo_);
      tree_->Branch("wt_l0_hi", &wt_l0_hi_);
      tree_->Branch("wt_l0_lo", &wt_l0_lo_);
      tree_->Branch("wt_trg_l0_hi", &wt_trg_l0_hi_);
      tree_->Branch("wt_trg_l0_lo", &wt_trg_l0_lo_);
      tree_->Branch("wt_p0_hi", &wt_p0_hi_);
      tree_->Branch("wt_p0_lo", &wt_p0_lo_);
      tree_->Branch("wt_p0_e_fake_hi", &wt_p0_e_fake_hi_);
      tree_->Branch("wt_p0_e_fake_lo", &wt_p0_e_fake_lo_);
      tree_->Branch("wt_p0_fake_err", &wt_p0_fake_err_);
      tree_->Branch("wt_p0_fake_bin", &wt_p0_fake_bin_);

      if (do_wg_gen_vars_) {
        tree_->Branch("is_wg_gen", &is_wg_gen_);
        tree_->Branch("gen_nparts", &gen_nparts_);

        tree_->Branch("gen_pdgid", &gen_pdgid_);

        tree_->Branch("gen_l0_match", &gen_l0_match_);
        tree_->Branch("gen_p0_match", &gen_p0_match_);

        tree_->Branch("lhe_l0_pt", &lhe_l0_pt_);
        tree_->Branch("lhe_l0_eta", &lhe_l0_eta_);
        tree_->Branch("lhe_p0_pt", &lhe_p0_pt_);
        tree_->Branch("lhe_p0_eta", &lhe_p0_eta_);
        // tree_->Branch("lhe_l0p0_dr", &lhe_l0p0_dr_);
        // tree_->Branch("lhe_p0j_dr", &lhe_p0j_dr_);
        // tree_->Branch("lhe_j_pt", &lhe_j_pt_);
      }
    }

    if (var_set_ >= 2) {
      tree_->Branch("run", &run_);

      tree_->Branch("l0_trg_2", &l0_trg_2_);

      tree_->Branch("p0_neiso", &p0_neiso_);
      tree_->Branch("p0_phiso", &p0_phiso_);
      tree_->Branch("p0_hovere", &p0_hovere_);
      tree_->Branch("p0_tight", &p0_tight_);
      tree_->Branch("p0_fsr", &p0_fsr_);
    }

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
  fns_["pileup_ratio_hi"] = std::shared_ptr<RooFunctor>(
    ws_->function("pileup_ratio_hi")->functor(ws_->argSet("pu_int")));
  fns_["pileup_ratio_lo"] = std::shared_ptr<RooFunctor>(
    ws_->function("pileup_ratio_lo")->functor(ws_->argSet("pu_int")));
  fns_["m_idisotrk_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_idisotrk_ratio")->functor(ws_->argSet("m_pt,m_eta")));
  fns_["m_idisotrk_ratio_err"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_idisotrk_ratio_err")->functor(ws_->argSet("m_pt,m_eta")));
  fns_["m_trg_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_trg_ratio")->functor(ws_->argSet("m_pt,m_eta")));
  fns_["m_trg_ratio_err"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_trg_ratio_err")->functor(ws_->argSet("m_pt,m_eta")));
  fns_["m_fake_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("m_fake_ratio")->functor(ws_->argSet("m_pt,m_eta")));

  fns_["e_gsfidiso_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_gsfidiso_ratio")->functor(ws_->argSet("e_pt,e_eta")));
  fns_["e_gsfidiso_ratio_err"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_gsfidiso_ratio_err")->functor(ws_->argSet("e_pt,e_eta")));
  fns_["e_trg_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_trg_ratio")->functor(ws_->argSet("e_pt,e_eta")));
  fns_["e_fake_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_fake_ratio")->functor(ws_->argSet("e_pt,e_eta")));

  fns_["p_id_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_id_ratio")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_psv_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_psv_ratio")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_id_ratio_err"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_id_ratio_err")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["e_p_fake_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_p_fake_ratio")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["e_p_fake_ratio_err"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_p_fake_ratio_err")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_fake_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_fake_ratio")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_fake_index"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_fake_index")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_fake_ratio_err"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_fake_ratio_err")->functor(ws_->argSet("p_pt,p_eta")));
  fns_["p_highpt_fake_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("p_highpt_fake_ratio")->functor(ws_->argSet("p_pt,p_eta")));

  if (rc_file_ != "") {
    rc_.init(rc_file_);
  }
  return 0;
  }

  int WGDataAnalysis::Execute(TreeEvent* event) {

    SetDefaults();

    if (check_is_wwg_) {
      unsigned n_neutrinos = 0;
      auto gen_parts = event->GetPtrVec<ac::GenParticle>("genParticles");
      for (auto const& p : gen_parts) {
        if (p->statusFlags().isPrompt() && p->status() == 1 && IsNeutrino(*p)) {
          ++n_neutrinos;
        }
      }
      if (n_neutrinos == 2) {
        return 1;
      }
    }

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    auto muons = event->GetPtrVec<ac::Muon>("muons");
    auto electrons = event->GetPtrVec<ac::Electron>("electrons");
    auto photons = event->GetPtrVec<ac::Photon>("photons");
    auto jets = event->GetPtrVec<ac::PFJet>("pfJets");

    double rho = event->Get<double>("fixedGridRhoFastjetAll");
    for (Photon *p : photons) {
      PhotonIsoCorrector(p, rho);
    }
    if (correct_p_energy_ >= 0) {
      if (year_ == 2018) {
        double corr_fac_2018 = correct_p_energy_ == 1 ? 1.002 : 0.998;
        for (Photon *p : photons) {
          p->setVector(p->vector() * corr_fac_2018);
        }
      } else {
        for (Photon *p : photons) {
          p->setVector(p->vector() * (p->energyCorrections().at(correct_p_energy_) / p->energy()));
        }
      }
    }
    if (correct_e_energy_ >= 0) {
      for (Electron *e : electrons) {
        e->setVector(e->vector() * (e->energyCorrections().at(correct_e_energy_) / e->energy()));
      }
    }
    if (correct_m_energy_ >= 0) {
      auto muons_to_correct = ac::copy_keep_if(muons, [](ac::Muon const* m) {
        return m->pt() > 15. && fabs(m->eta()) < 2.4;
      });
      if (is_data_) {
        for (Muon *m : muons_to_correct) {
          m->setVector(m->vector() * rc_.kScaleDT(m->charge(), m->pt(), m->eta(), m->phi()));
        }
      } else if (muons_to_correct.size() > 0) {
        auto gen_parts = event->GetPtrVec<ac::GenParticle>("genParticles");
        auto gen_muons = ac::copy_keep_if(gen_parts, [&](ac::GenParticle *p) {
          return ac::IsMuon(*p) && p->status() == 1 && p->pt() > 10.0;
        });
        for (Muon *m : muons_to_correct) {
          auto matched_muons = ac::copy_keep_if(gen_muons, [&](ac::GenParticle *p) {
            return ac::DeltaR(m, p) < 0.3;
          });
          boost::range::sort(matched_muons, DescendingPt);
          if (matched_muons.size() > 0) {
            double sf =
                rc_.kSpreadMC(m->charge(), m->pt(), m->eta(), m->phi(), matched_muons.at(0)->pt());
            if (correct_m_energy_ >= 1) {
              // Do systematic variations
              double delta_sf = rc_.kSpreadMCerror(m->charge(), m->pt(), m->eta(), m->phi(), matched_muons.at(0)->pt());
              if (correct_m_energy_ == 1) sf += delta_sf;
              if (correct_m_energy_ == 2) sf -= delta_sf;
            }
            m->setVector(m->vector() * sf);
          }
        }
      }

    }
    auto mets = event->GetPtrVec<ac::Met>("pfType1Met");
    auto puppi_mets = event->GetPtrVec<ac::Met>("puppiMet");
    // auto tk_mets = event->GetPtrVec<ac::Met>("trackMet");

    // Sub-process classification - maybe move this into a separate module
    /*
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
    */

    auto pre_muons = ac::copy_keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && MuonPFIso(m) < 0.15 &&
             fabs(m->dxy()) < 0.05 && fabs(m->dz()) < 0.2;
    });

    auto pre_muons_loose_iso = ac::copy_keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() &&
             fabs(m->dxy()) < 0.05 && fabs(m->dz()) < 0.2;
    });

    auto pre_electrons = ac::copy_keep_if(electrons, [](ac::Electron const* e) {
      return e->pt() > 35. && fabs(e->scEta()) < 2.5 && e->isCutBasedMediumElectron() &&
             ElectronIsoFall17V2(e, 2) && (fabs(e->scEta()) < 1.4442 || fabs(e->scEta()) > 1.566) &&
             ElectronIPCuts(e);
    });

    auto pre_electrons_loose_iso = ac::copy_keep_if(electrons, [](ac::Electron const* e) {
      return e->pt() > 35. && fabs(e->scEta()) < 2.5 && e->isCutBasedMediumElectron() &&
             (fabs(e->scEta()) < 1.4442 || fabs(e->scEta()) > 1.566) && ElectronIPCuts(e);
    });
    // At this stage apply the medium Photon ID without the charged iso cut
    auto pre_photons = ac::copy_keep_if(photons, [&](ac::Photon const* p) {
      return p->pt() > 30. && fabs(p->scEta()) < 2.5 &&
             (fabs(p->scEta()) < 1.4442 || fabs(p->scEta()) > 1.566) && PhotonIDIso(p, year_, 1, false, false);
    });

    bool super_loose = true;
    double super_loose_threshold = 100.;
    auto super_loose_photons = ac::copy_keep_if(photons, [&](ac::Photon const* p) {
      return p->pt() > 100. && fabs(p->scEta()) < 2.5 &&
             (fabs(p->scEta()) < 1.4442 || fabs(p->scEta()) > 1.566) && super_loose && p->pt() > super_loose_threshold && p->hadTowOverEm() < 0.15;
    });

    // TODO: update with proper loose ID
    auto veto_muons = ac::copy_keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 10. && fabs(m->eta()) < 2.4 && fabs(m->dxy()) < 0.05 &&
             fabs(m->dz()) < 0.2;
    });

    auto veto_electrons = ac::copy_keep_if(electrons, [](ac::Electron const* e) {
      return e->pt() > 10. && fabs(e->scEta()) < 2.5 && e->isCutBasedLooseElectron() &&
             ElectronIsoFall17V2(e, 1) && (fabs(e->scEta()) < 1.4442 || fabs(e->scEta()) > 1.566) &&
             ElectronIPCuts(e);
    });

    auto alt_veto_electrons = veto_electrons;
    auto alt_veto_muons = veto_muons;

    auto qcd_jets = ac::copy_keep_if(jets, [](ac::PFJet const* j) {
      return j->pt() > 30. && fabs(j->eta()) < 2.5 && j->passesJetID();
    });

    boost::range::sort(pre_muons, DescendingPt);
    boost::range::sort(pre_electrons, DescendingPt);
    boost::range::sort(pre_muons_loose_iso, DescendingPt);
    boost::range::sort(pre_electrons_loose_iso, DescendingPt);
    boost::range::sort(pre_photons, DescendingPt);
    boost::range::sort(super_loose_photons, DescendingPt);
    boost::range::sort(qcd_jets, DescendingPt);

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
    if (l0) l0_nominal_ = true;

    // The event doesn't have a viable lepton, we can try and fallback
    // to a loose iso one
    if (!l0) {
      if (pre_muons_loose_iso.size() > 0) {
        m0 = pre_muons_loose_iso[0];
        pre_muons.push_back(m0);
      }
      if (pre_electrons_loose_iso.size() > 0) {
        e0 = pre_electrons_loose_iso[0];
        pre_electrons.push_back(e0);
      }
      // And apply the preference logic again
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
    }

    if (m0) {
      l0_pdgid_ = 13;
    } else if (e0) {
      l0_pdgid_ = 11;
    }


    if (l0) {
      ac::keep_if(pre_photons, [&](ac::Photon const* p) {
        return DeltaR(p, l0) > 0.7;
      });
      ac::keep_if(super_loose_photons, [&](ac::Photon const* p) {
        return DeltaR(p, l0) > 0.7;
      });
      ac::keep_if(qcd_jets, [&](ac::PFJet const* j) {
        return DeltaR(j, l0) > 0.7;
      });
      ac::keep_if(veto_electrons, [&](ac::Electron const* e) {
        return DeltaR(e, l0) > 0.3;
      });
      ac::keep_if(veto_muons, [&](ac::Muon const* m) {
        return DeltaR(m, l0) > 0.3;
      });
      ac::keep_if(alt_veto_electrons, [&](ac::Electron const* e) {
        return e != l0;
      });
      ac::keep_if(alt_veto_muons, [&](ac::Muon const* m) {
        return m != l0;
      });
    }

    if (pre_photons.size() == 0) {
      pre_photons = super_loose_photons;
    }

    ac::Photon* p0 = pre_photons.size() ? pre_photons[0] : nullptr;

    ac::PFJet* j0 = qcd_jets.size() ? qcd_jets[0] : nullptr;

    unsigned use_met = 0;
    if (shift_met_ >= 0) {
      use_met = shift_met_;
    }

    ac::Met* met = mets.at(use_met);
    ac::Met* puppi_met = puppi_mets.at(use_met);

    // Here are the conditions for skipping this event and not writing anything to the tree
    bool wg_presel = (l0 != nullptr) && (p0 != nullptr); // any
    // bool lj_presel = (l0 != nullptr) && (j0 != nullptr) && puppi_met->pt() < 30.0 && MT(l0, puppi_met) < 30.0;
    bool lj_presel = (l0 != nullptr) && (j0 != nullptr);
    bool mm_presel = (pre_muons.size() >= 2 && l0_nominal_);
    bool ee_presel = (pre_electrons.size() >= 2 && l0_nominal_);
    bool reco_event = only_wg_ ? wg_presel : (wg_presel || mm_presel || ee_presel || lj_presel);
    if (do_presel_ && !reco_event) {
      return 1;
    }

    n_pre_m_ = pre_muons.size();
    n_pre_e_ = pre_electrons.size();
    n_pre_p_ = pre_photons.size();

    n_veto_m_ = veto_muons.size();
    n_veto_e_ = veto_electrons.size();
    n_alt_veto_m_ = alt_veto_muons.size();
    n_alt_veto_e_ = alt_veto_electrons.size();

    n_qcd_j_ = qcd_jets.size();

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
        auto lhe_partons = ac::copy_keep_if(lhe_parts, [](ac::GenParticle *p) {
          unsigned apdgid = std::abs(p->pdgId());
          return p->status() == 1 && ((apdgid >= 1 && apdgid <= 6) || apdgid == 21);
        });
        boost::range::sort(lhe_partons, [&](ac::GenParticle* x1, ac::GenParticle* x2) {
          return DeltaR(x1, parts.lhe_pho) < DeltaR(x2, parts.lhe_pho);
        });
        if (lhe_partons.size() > 0) {
          lhe_p0j_dr_ = DeltaR(parts.lhe_pho, lhe_partons[0]);
          lhe_j_pt_ = lhe_partons[0]->pt();
        }
        double frixione_sum = 0.;
        double frixione_dr = 0.4;
        for (auto const& ip : lhe_partons) {
          double dr = DeltaR(ip, parts.lhe_pho);
          if (dr >= frixione_dr) {
            break;
          }
          frixione_sum += ip->pt();
          if (frixione_sum > (parts.lhe_pho->pt() * ((1. - std::cos(dr)) / (1. - std::cos(frixione_dr))))) {
            lhe_frixione_ = false;
          }
        }
      }
      if (parts.ok) {
        is_wg_gen_ = true;
        WGSystem gen_sys = ProduceWGSystem(*parts.gen_lep, *gen_met, *parts.gen_pho, true, rng, false);
        WGSystem gen_true_sys = ProduceWGSystem(*parts.gen_lep, *parts.gen_neu, *parts.gen_pho, false, rng, false);
        gen_p0_pt_ = parts.gen_pho->pt();
        gen_p0_eta_ = parts.gen_pho->eta();
        gen_l0_q_ = parts.gen_lep->charge();
        gen_phi_ = gen_sys.Phi(parts.gen_lep->charge() > 0);
        gen_phi_f_ = gen_sys.SymPhi(parts.gen_lep->charge() > 0);
        gen_l0_pt_ = parts.gen_lep->pt();
        gen_l0_eta_ = parts.gen_lep->eta();
        gen_met_ = gen_met->pt();
        gen_l0p0_dr_ = ac::DeltaR(parts.gen_lep, parts.gen_pho);
        true_phi_ = gen_true_sys.Phi(parts.gen_lep->charge() > 0);
        true_phi_f_ = gen_true_sys.SymPhi(parts.gen_lep->charge() > 0);
        gen_wg_M_ = (parts.gen_lep->vector() + parts.gen_neu->vector() + parts.gen_pho->vector()).M();
        gen_mt_cluster_ = ac::MTcluster(parts.gen_lep, parts.gen_pho, gen_met);
        // Now try and match to the reco objects, if they exist
        if (l0 && ac::DeltaR(l0, parts.gen_lep) < 0.3) {
          gen_l0_match_ = true;
        }
        if (p0 && ac::DeltaR(p0, parts.gen_pho) < 0.3) {
          gen_p0_match_ = true;
        }
      }

      // const auto rivet = event->GetPtr<WGammaRivetVariables>("rivetVariables");
      // bool fid = is_wg_gen_ && gen_l0_pt_ > 30. && std::abs(gen_l0_eta_) < 2.5 && gen_p0_pt_ > 30. && std::abs(gen_p0_eta_) < 2.5 && gen_l0p0_dr_ > 0.7 && lhe_frixione_ && gen_met_ > 40.;
      // fid = fid && lhe_frixione_;
      // bool rfid = rivet->is_wg_gen && rivet->l0_pt > 30. && std::abs(rivet->l0_eta) < 2.5 && rivet->p0_pt > 30. && std::abs(rivet->p0_eta) < 2.5 && rivet->l0p0_dr > 0.7 && rivet->p0_frixione && rivet->met_pt > 40.;
      // if (fid != rfid) {
      //   std::cout << "is_wg_gen: " << is_wg_gen_ << "\t" << rivet->is_wg_gen << "\n";
      //   std::cout << "l0_pt: " << gen_l0_pt_ << "\t" << rivet->l0_pt << "\n";
      //   std::cout << "l0_eta: " << gen_l0_eta_ << "\t" << rivet->l0_eta << "\n";
      //   std::cout << "l0_q: " << gen_l0_q_ << "\t" << rivet->l0_q << "\n";
      //   std::cout << "l0_pdgid: " << gen_pdgid_ << "\t" << rivet->l0_abs_pdgid << "\n";
      //   std::cout << "p0_pt: " << gen_p0_pt_ << "\t" << rivet->p0_pt << "\n";
      //   std::cout << "l0p0_dr: " << gen_l0p0_dr_ << "\t" << rivet->l0p0_dr << "\n";
      //   std::cout << "p0_frixione: " << lhe_frixione_ << "\t" << rivet->p0_frixione << "\n";
      //   std::cout << "met_pt: " << gen_met_ << "\t" << rivet->met_pt << "\n";
      // }
      // if (fid && rfid) {
      //   std::cout << "mt_cluster: " << gen_mt_cluster_ << "\t" << rivet->mt_cluster << "\n";
      // }
    }

    if (reco_event) {
      run_ = info->run();
      n_vtx_ = info->numVertices();
      metfilters_ = info->metfilters().to_ulong();

      met_ = met->pt();
      met_phi_ = met->phi();
      puppi_met_ = puppi_met->pt();
      puppi_met_phi_ = puppi_met->phi();
    }

    if (l0) {
      l0_pt_ = l0->pt();
      l0_eta_ = l0->eta();
      l0_phi_ = l0->phi();
      l0_q_ = l0->charge();
      l0met_mt_ = ac::MT(l0, puppi_met);

      if (m0) {
        l0_iso_ = MuonPFIso(m0);
        l0_tkiso_ = m0->pfIsoSumChargedHadronPt() / m0->pt();
        l0_medium_ = m0->isMediumMuon();
        l0_tight_ = m0->isTightMuon();
      } else if (e0) {
        l0_iso_ = e0->relativeEAIso();
        l0_medium_ = e0->isCutBasedMediumElectron();
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
              IsFilterMatchedDR(m0, trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3) ||
              IsFilterMatchedDR(m0, trg_objs_tk, filters_IsoTkMu24_.Lookup(trg_lookup), 0.3);
        } else if (year_ == 2017) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu27");
          l0_trg_ = IsFilterMatchedDR(m0, trg_objs, filters_IsoMu27_.Lookup(trg_lookup), 0.3);
        } else if (year_ == 2018) {
          auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
          l0_trg_ = IsFilterMatchedDR(m0, trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3);
        }
        auto const& trg_objs_Mu50 = event->GetPtrVec<TriggerObject>("triggerObjects_Mu50");
        l0_trg_2_ = IsFilterMatchedDR(m0, trg_objs_Mu50, filters_Mu50_.Lookup(trg_lookup), 0.3);
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

    if (m0 && pre_muons.size() >= 2) {
      ac::Muon* m1 = pre_muons[1];
      l1_pt_ = m1->pt();
      l1_eta_ = m1->eta();
      l1_phi_ = m1->phi();
      l1_iso_ = MuonPFIso(m1);

      l0l1_M_ = (m0->vector() + m1->vector()).M();
      l0l1_dr_ = DeltaR(m0, m1);
      l0l1_os_ = m0->charge() != m1->charge();
      l0l1_pt_ = (m0->vector() + m1->vector()).pt();
    } else if (e0 && pre_electrons.size() >= 2) {
      ac::Electron* e1 = pre_electrons[1];
      l1_pt_ = e1->pt();
      l1_eta_ = e1->eta();
      l1_phi_ = e1->phi();
      l1_iso_ = e1->relativeEAIso();

      l0l1_M_ = (e0->vector() + e1->vector()).M();
      l0l1_dr_ = DeltaR(e0, e1);
      l0l1_os_ = e0->charge() != e1->charge();
      l0l1_pt_ = (e0->vector() + e1->vector()).pt();
    }

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
      p0_loose_ = p0->isLooseIdPhoton();
      p0_medium_ = p0->isMediumIdPhoton();
      p0_tight_ = p0->isTightIdPhoton();

      // if (n_vm_ >= 1) {
      //   boost::range::sort(veto_muons, [&](ac::Muon* x1, ac::Muon* x2) {
      //     return DeltaR(x1, p0) < DeltaR(x2, p0);
      //   });
      //   vm_p0_dr_ = ac::DeltaR(veto_muons[0], p0);
      // }
    }

    if (j0) {
      j0_pt_ = j0->pt();
      j0_eta_ = j0->eta();
      if (l0) {
        l0j0_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(l0->vector(), j0->vector());
      }
    }

    if (p0 && l0) {
      l0p0_dr_ =  ac::DeltaR(l0, p0);
      l0p0_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(l0->vector(), p0->vector());
      l0p0_M_ = (l0->vector() + p0->vector()).M();

      // special transverse mass
      mt_cluster_ = ac::MTcluster(l0, p0, met);

      WGSystem reco_sys = ProduceWGSystem(*l0, *met, *p0, true, rng, false);
      reco_phi_ = reco_sys.Phi(l0->charge());
      reco_phi_f_ = reco_sys.SymPhi(l0->charge());

      // WGSystem reco_tk_sys = ProduceWGSystem(*l0, *tk_met, *p0, true, rng, false);
      // reco_tk_phi_ = reco_tk_sys.Phi(l0->charge());
      // reco_tk_phi_f_ = reco_tk_sys.SymPhi(l0->charge());

      WGSystem reco_puppi_sys = ProduceWGSystem(*l0, *puppi_met, *p0, true, rng, false);
      reco_puppi_phi_ = reco_puppi_sys.Phi(l0->charge());
      reco_puppi_phi_f_ = reco_puppi_sys.SymPhi(l0->charge());

      wt_p0_fake_ = RooFunc(fns_["p_fake_ratio"], {p0_pt_, p0->scEta()});
      if (super_loose && p0_pt_ > super_loose_threshold) {
        wt_p0_highpt_fake_ = RooFunc(fns_["p_highpt_fake_ratio"], {p0_pt_, p0->scEta()});
      } else {
        wt_p0_highpt_fake_ = 0.;
      }
      wt_p0_fake_err_ = RooFunc(fns_["p_fake_ratio_err"], {p0_pt_, p0->scEta()});
      wt_p0_fake_bin_ = RooFunc(fns_["p_fake_index"], {p0_pt_, p0->scEta()});
      // wt_p0_fake_lo_ = 1.0 - RooFunc(fns_["p_fake_ratio_err"], {p0_pt_, p0->scEta()});
      if (e0) {
        wt_l0_fake_ = RooFunc(fns_["e_fake_ratio"], {l0_pt_, l0_eta_});
      }
      if (m0) {
        wt_l0_fake_ = RooFunc(fns_["m_fake_ratio"], {l0_pt_, l0_eta_});
      }
    }

    if (!is_data_) {
      // wt_def_ = info->nominalGenWeight() >= 0. ? +1. : -1.; // temporary until new ntuples fix EventInfo bug
      wt_def_ = info->totalWeight();
      auto const& pu_info = event->GetPtrVec<PileupInfo>("pileupInfo");
      for (PileupInfo const* pu : pu_info) {
        if (pu->bunchCrossing() == 0) {
          wt_pu_ = RooFunc(fns_["pileup_ratio"], {pu->trueNumInteractions()});
          wt_pu_hi_ = RooFunc(fns_["pileup_ratio_hi"], {pu->trueNumInteractions()}) / wt_pu_;
          wt_pu_lo_ = RooFunc(fns_["pileup_ratio_lo"], {pu->trueNumInteractions()}) / wt_pu_;
          break;
        }
      }
      if (year_ == 2016 || year_ == 2017) {
        wt_pf_ = event->Get<double>("NonPrefiringProb");
        wt_pf_lo_ = event->Get<double>("NonPrefiringProbUp") / wt_pf_;
        wt_pf_hi_ = event->Get<double>("NonPrefiringProbDown") / wt_pf_;
      }
      if (m0) {
        wt_l0_ = RooFunc(fns_["m_idisotrk_ratio"], {l0_pt_, l0_eta_});
        wt_l0_hi_ = 1.0 + RooFunc(fns_["m_idisotrk_ratio_err"], {l0_pt_, l0_eta_});
        wt_l0_lo_ = 1.0 - RooFunc(fns_["m_idisotrk_ratio_err"], {l0_pt_, l0_eta_});
        wt_trg_l0_ = RooFunc(fns_["m_trg_ratio"], {l0_pt_, l0_eta_});
        wt_trg_l0_hi_ = 1.0 + RooFunc(fns_["m_trg_ratio_err"], {l0_pt_, l0_eta_});
        wt_trg_l0_lo_ = 1.0 - RooFunc(fns_["m_trg_ratio_err"], {l0_pt_, l0_eta_});
      }
      if (e0) {
        // For electrons the SFs are binned in scEta instead of eta
        wt_l0_ = RooFunc(fns_["e_gsfidiso_ratio"], {l0_pt_, e0->scEta()});
        wt_l0_hi_ = 1.0 + RooFunc(fns_["e_gsfidiso_ratio_err"], {l0_pt_, e0->scEta()});
        wt_l0_lo_ = 1.0 - RooFunc(fns_["e_gsfidiso_ratio_err"], {l0_pt_, e0->scEta()});
        wt_trg_l0_ = RooFunc(fns_["e_trg_ratio"], {l0_pt_, e0->scEta()});
      }
      if (m0 && pre_muons.size() >= 2) {
        wt_l1_ = RooFunc(fns_["m_idisotrk_ratio"], {l1_pt_, l1_eta_});
      }
      if (e0 && pre_electrons.size() >= 2) {
        wt_l1_ = RooFunc(fns_["e_gsfidiso_ratio"], {l1_pt_, pre_electrons[1]->scEta()});
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

          if (check_is_wwg_) {
            int p0_idx = prompt_gen_photons[0]->index();
            for (auto const& p : genparts) {
              if (ac::contains(p->daughters(), p0_idx)) {
                // These flags should only select photons that come from V->qq(q->qgamma)
                if (p->status() == 23 && std::abs(p->pdgId()) >= 1 && std::abs(p->pdgId()) <= 6 && p->statusFlags().isLastCopyBeforeFSR()) {
                  gen_proc_ = 1;
                }
              }
            }
          }
        } else {
          // So not any of the above, is it at least matched to a genJet with pT > 15 GeV?
          auto genjets = event->GetPtrVec<ac::Candidate>("genJets");
          auto genjet_dr = ac::copy_keep_if(genjets, [&](ac::Candidate *p) {
            return p->pt() > 15.0 && DeltaR(p, p0) < 0.5;
          });
          if (genjet_dr.size() == 0) {
            p0_truth_ = 7; // Pileup photon
          }

        }
        // if (p0_truth_ == 0 && p0_medium_) {
        //   std::cout << ">>>> Fake hadronic photon found: " << p0->vector() << "\n";
        // }
        // if (p0_truth_ == 0 && p0_medium_) {
        //   auto gen_dr = ac::copy_keep_if(genparts, [&](ac::GenParticle *p) {
        //     return DeltaR(p, p0) < 0.5;
        //   });

        //   std::cout << ">>>> Fake photon found: " << p0->vector() << "\n";
        //   std::cout << ">>>> Nearby GenParticles:\n";
        //   for (auto const& gpart : gen_dr) {
        //     gpart->Print();
        //   }
        //   std::cout << ">>>> Nearby GenJets:\n";
        //   for (auto const& gpart : genjet_dr) {
        //     gpart->Print();
        //   }

        // }
        // for (auto const& p_pho : prompt_gen_photons) {
        //   for (auto const& p : genparts) {
        //     if (ac::IsChargedLepton(*p) && ac::contains(p->daughters(), p_pho->index())) {
        //       p0_fsr_ = true;
        //       break;
        //     }
        //   }
        //   if (p0_fsr_) break;
        // }
        if (ac::contains({1, 2, 4, 5}, p0_truth_)) {
          wt_p0_ = RooFunc(fns_["p_id_ratio"], {p0_pt_, p0->scEta()});
          wt_p0_hi_ = 1.0 + RooFunc(fns_["p_id_ratio_err"], {p0_pt_, p0->scEta()});
          wt_p0_lo_ = 1.0 - RooFunc(fns_["p_id_ratio_err"], {p0_pt_, p0->scEta()});
        }
        if (ac::contains({1, 4, 5}, p0_truth_) && !p0->hasPixelSeed() && p0->passElectronVeto()) {
          wt_p0_ *= RooFunc(fns_["p_psv_ratio"], {p0_pt_, p0->scEta()});
        }
        if (ac::contains({2}, p0_truth_) && !p0->hasPixelSeed() && p0->passElectronVeto()) {
          wt_p0_e_fake_ = RooFunc(fns_["e_p_fake_ratio"], {p0_pt_, p0->eta()});
          wt_p0_e_fake_hi_ = 1.0 + RooFunc(fns_["e_p_fake_ratio_err"], {p0_pt_, p0->eta()});
          wt_p0_e_fake_lo_ = 1.0 - RooFunc(fns_["e_p_fake_ratio_err"], {p0_pt_, p0->eta()});
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

      if (check_gen_mll_) {
        ROOT::Math::PtEtaPhiMVector dy_parts;
        unsigned n_leps = 0;
        for (auto const& p : genparts) {
          if (p->statusFlags().isPrompt() && IsChargedLepton(*p) && p->statusFlags().isLastCopy()) {
            dy_parts += p->vector();
            ++n_leps;
          }
        }
        if (n_leps >= 2) {
          gen_mll_ = dy_parts.M();
        }
      }

      if (scale_weights_ > 0) {
        auto const& scale_weights = event->Get<std::vector<double>>("scale_weights");
        wt_sc_0_ = scale_weights.at(0);
        wt_sc_1_ = scale_weights.at(1);
        wt_sc_2_ = scale_weights.at(2);
        wt_sc_3_ = scale_weights.at(3);
        wt_sc_4_ = scale_weights.at(4);
        wt_sc_5_ = scale_weights.at(5);
      }

      int ipdf = 0;
      for (unsigned ipdfset = 0; ipdfset < pdf_begin_.size(); ++ipdfset) {
        int pdf_begin = pdf_begin_.at(ipdfset);
        int pdf_end = pdf_end_.at(ipdfset);
        int npdf = pdf_end - pdf_begin + 1;
        for (int ip = 0; ip < npdf; ++ip) {
          if (info->lheWeights().count(pdf_begin + ip)) {
            wt_pdf_[ipdf] = 1. + 2. * info->lheWeights().at(pdf_begin + ip);
          } else {
            wt_pdf_[ipdf] = 0.;
          }
          ++ipdf;
          // std::cout << " - pdf " << ipdf << "/" << pdf_begin_ + ipdf << "\t" << wt_pdf_[ipdf] << "\n";
        }
      }

      if (ps_weights_ > 0) {
        auto const& gen_weights = info->genWeights();
        if (gen_weights.size() >= 14) {
          wt_isr_hi_ = gen_weights[6];
          wt_isr_lo_ = gen_weights[8];
          wt_fsr_hi_ = gen_weights[7];
          wt_fsr_lo_ = gen_weights[9];
        } else {
          // std::cout << "Not enough elements!\n";
        }
      }
    }
    CompressVars();
    //if (nproc == 10000) {
    //  tree_->OptimizeBaskets(10E6, 1.1, "d");
    //}
    tree_->Fill();
    //++nproc;
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
    n_alt_veto_m_ = 0;
    n_alt_veto_e_ = 0;
    n_qcd_j_ = 0;
    l0_pt_ = 0.;
    l0_eta_ = 0.;
    l0_phi_ = 0.;
    l0_iso_ = 0.;
    l0_nominal_ = false;
    l0_tkiso_ = 0.;
    l0_medium_ = false;
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
    l0l1_pt_ = 0.;
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
    p0_loose_ = false;
    p0_fsr_ = false;
    p0_tight_ = false;
    p0_truth_ = 0;
    j0_pt_ = 0.;
    j0_eta_ = 0.;
    l0j0_dphi_ = 0.;
    met_ = 0.;
    met_phi_ = 0.;
    tk_met_ = 0.;
    tk_met_phi_ = 0.;
    puppi_met_ = 0.;
    puppi_met_phi_ = 0.;
    l0met_mt_ = 0.;
    l0p0_dr_ = 0.;
    l0p0_dphi_ = 0.;
    l0p0_M_ = 0.;
    mt_cluster_ = 0.;
    reco_phi_ = 0.;
    reco_phi_f_ = 0.;
    reco_tk_phi_ = 0.;
    reco_tk_phi_f_ = 0.;
    reco_puppi_phi_ = 0.;
    reco_puppi_phi_f_ = 0.;
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
    wt_l0_fake_ = 1.;
    wt_p0_e_fake_ = 1.;
    wt_sc_0_ = 1.;
    wt_sc_1_ = 1.;
    wt_sc_2_ = 1.;
    wt_sc_3_ = 1.;
    wt_sc_4_ = 1.;
    wt_sc_5_ = 1.;
    wt_isr_hi_ = 1.;
    wt_isr_lo_ = 1.;
    wt_fsr_hi_ = 1.;
    wt_fsr_lo_ = 1.;
    wt_pu_hi_ = 1.;
    wt_pu_lo_ = 1.;
    wt_pf_hi_ = 1.;
    wt_pf_lo_ = 1.;
    wt_l0_hi_ = 1.;
    wt_l0_lo_ = 1.;
    wt_trg_l0_hi_ = 1.;
    wt_trg_l0_lo_ = 1.;
    wt_p0_hi_ = 1.;
    wt_p0_lo_ = 1.;
    wt_p0_e_fake_hi_ = 1.;
    wt_p0_e_fake_lo_ = 1.;
    wt_p0_fake_err_ = 1.;
    wt_p0_fake_bin_ = 0;
    is_wg_gen_ = false;
    gen_pdgid_ = 0;
    gen_nparts_ = 0;
    gen_l0_match_ = false;
    gen_p0_match_ = false;
    lhe_l0_pt_ = 0.;
    lhe_l0_eta_ = 0.;
    lhe_p0_pt_ = 0.;
    lhe_p0_eta_ = 0.;
    lhe_l0p0_dr_ = 0.;
    lhe_p0j_dr_ = -1.;
    lhe_j_pt_ = -1.;
    lhe_frixione_ = true;
    gen_p0_pt_ = 0.;
    gen_p0_eta_ = 0.;
    gen_phi_ = 0.;
    gen_phi_f_ = 0.;
    gen_l0_q_ = 0;
    gen_l0_pt_ = 0.;
    gen_wg_M_ = 0.;
    gen_mt_cluster_ = 0.;
    gen_l0_eta_ = 0.;
    gen_met_ = 0.;
    gen_l0p0_dr_ = 0.;
    true_phi_ = 0.;
    true_phi_f_ = 0.;
    gen_mll_ = -1.0;
    // gen_n2_pt_ = 0.;
  }

  void WGDataAnalysis::CompressVars() {
    for (float* var :
         {&l0_pt_,     &l0_iso_,    &l0met_mt_,  &l1_pt_,     &l1_iso_,      &met_,
          &puppi_met_, &tk_met_,    &l0met_mt_,  &l0l1_M_,    &l0l1_pt_,   &l0l1_dr_,     &lhe_l0_pt_,
          &lhe_p0_pt_, &gen_p0_pt_, &gen_l0_pt_, &gen_wg_M_,  &gen_met_,   &gen_l0p0_dr_, &p0_pt_,
          &j0_pt_, &mt_cluster_, &gen_mt_cluster_, &gen_mll_,
          &p0_chiso_,  &p0_neiso_,  &p0_phiso_,  &p0_hovere_, &p0_sigma_,    &l0p0_dr_,
          &l0p0_dphi_, &l0p0_M_, &wt_sc_0_, &wt_sc_1_, &wt_sc_2_, &wt_sc_3_, &wt_sc_4_, &wt_sc_5_,
          &wt_isr_lo_, &wt_isr_hi_, &wt_fsr_lo_, &wt_fsr_hi_}) {
      *var = reduceMantissaToNbitsRounding(*var, 10);
    }
    for (float* var :
         {&l0_eta_,          &l0_phi_,     &l1_eta_,     &l1_phi_,     &met_phi_,
          &puppi_met_phi_,   &tk_met_phi_, &lhe_l0_eta_, &lhe_p0_eta_, &gen_p0_eta_, &gen_phi_,
          &gen_phi_f_,       &gen_l0_eta_, &true_phi_,   &true_phi_f_, &p0_eta_,
          &j0_eta_, &l0j0_dphi_,
          &p0_phi_,          &l0p0_dphi_,  &reco_phi_,   &reco_phi_f_, &reco_puppi_phi_,
          &reco_puppi_phi_f_, &reco_tk_phi_, &reco_tk_phi_f_}) {
      *var = reduceMantissaToNbitsRounding(*var, 8);
    }
  }

  int WGDataAnalysis::PostAnalysis() {
    return 0;
  }

  void WGDataAnalysis::PrintInfo() {}

  void WGDataAnalysis::PhotonIsoCorrector(ac::Photon *p, double rho) {
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
