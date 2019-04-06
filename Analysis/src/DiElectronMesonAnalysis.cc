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
#include "Acorn/Analysis/interface/DiElectronMesonAnalysis.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"

namespace ac {

DiElectronMesonAnalysis::DiElectronMesonAnalysis(std::string const& name)
    : ModuleBase(name), fs_(nullptr), year_(2016), is_data_(true) {}

DiElectronMesonAnalysis::~DiElectronMesonAnalysis() { ; }

int DiElectronMesonAnalysis::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("DiElectronMesonAnalysis", "DiElectronMesonAnalysis");
    tree_->Branch("electron_overlaps",&electron_overlaps_);
    tree_->Branch("muon_overlaps",&muon_overlaps_);
    tree_->Branch("pt_1", &pt_1_);
    tree_->Branch("pt_2", &pt_2_);
    tree_->Branch("eta_1", &eta_1_);
    tree_->Branch("eta_2", &eta_2_);
    tree_->Branch("m_ll", &m_ll_);
    tree_->Branch("pt_ll", &pt_ll_);
    tree_->Branch("dr_ll", &dr_ll_);
    tree_->Branch("trg_1", &trg_1_);
    tree_->Branch("trg_2", &trg_2_);
    tree_->Branch("wt_pu", &wt_pu_);
    tree_->Branch("wt_1", &wt_1_);
    tree_->Branch("wt_2", &wt_2_);
    tree_->Branch("wt_trg1", &wt_trg1_);
    tree_->Branch("wt_trg2", &wt_trg2_);
    tree_->Branch("wt_rhoiso", &wt_rhoiso_);
    tree_->Branch("highestpt_pair_id_1", &highestpt_pair_id_1_);
    tree_->Branch("highestpt_pair_id_2", &highestpt_pair_id_2_);
    tree_->Branch("highestpt_pair_dR", &highestpt_pair_dR_);
    tree_->Branch("highestpt_pair_eta", &highestpt_pair_eta_);
    tree_->Branch("highestpt_pair_phi", &highestpt_pair_phi_);
    tree_->Branch("highestpt_pair_1_eta", &highestpt_pair_1_eta_);
    tree_->Branch("highestpt_pair_1_phi", &highestpt_pair_1_phi_);
    tree_->Branch("highestpt_pair_1_pt", &highestpt_pair_1_pt_);
    tree_->Branch("highestpt_pair_2_eta", &highestpt_pair_2_eta_);
    tree_->Branch("highestpt_pair_2_phi", &highestpt_pair_2_phi_);
    tree_->Branch("highestpt_pair_2_pt", &highestpt_pair_2_pt_);
    tree_->Branch("highestpt_pair_mass", &highestpt_pair_mass_);
    tree_->Branch("highestpt_pair_pt", &highestpt_pair_pt_);
    tree_->Branch("highestpt_pair_iso", &highestpt_pair_iso_);
    tree_->Branch("highestpt_pair_looser_iso", &highestpt_pair_looser_iso_);
    tree_->Branch("highestpt_pair_reco_higgs_mass", &highestpt_pair_reco_higgs_mass_);
    tree_->Branch("highestpt_pair_reco_higgs_pt", &highestpt_pair_reco_higgs_pt_);
    tree_->Branch("Zrho_dphi",&Zrho_dphi_);
    tree_->Print();
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
  fns_["e_gsfidiso_ratio"] = std::shared_ptr<RooFunctor>(
    ws_->function("e_gsfidiso_ratio")->functor(ws_->argSet("e_pt,e_eta")));
  fns_["rhoiso_ratio_etainc"] = std::shared_ptr<RooFunctor>(
    ws_->function("rhoiso_ratio_etainc")->functor(ws_->argSet("rho_pt,rho_eta")));

  return 0;
  }

  int DiElectronMesonAnalysis::Execute(TreeEvent* event) {

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    std::vector<ac::Muon *> muons = event->GetPtrVec<ac::Muon>("muons");
    std::vector<ac::Muon *> veto_muons = event->GetPtrVec<ac::Muon>("muons");
    std::vector<ac::Muon *> all_muons = event->GetPtrVec<ac::Muon>("muons");
    std::vector<ac::Electron *> electrons = event->GetPtrVec<ac::Electron>("electrons");
    std::vector<ac::Electron *> veto_electrons = event->GetPtrVec<ac::Electron>("electrons");
    std::vector<ac::Electron *> all_electrons = event->GetPtrVec<ac::Electron>("electrons");
    std::vector<ac::Track *> tracks = event->GetPtrVec<ac::Track>("Tracks");
    std::vector<ac::Track *> tracksforiso = event->GetPtrVec<ac::Track>("TracksForIso");
    std::vector<int> higgs_daughters;
    std::vector<int> rho_daughters;
    std::vector<GenParticle *> pions_from_meson;



    // Apply pT and ID cuts
    // The medium ID already includes iso?
    ac::keep_if(veto_muons, [](ac::Muon const* m) {
      return m->pt() > 5. && fabs(m->eta()) < 2.4 && m->isMediumMuon();
    });

    ac::keep_if(veto_electrons, [](ac::Electron const* e) {
      return e->pt() > 5. && fabs(e->eta()) < 2.1 && e->isMVAwp90Electron();
    });

    ac::keep_if(electrons, [](ac::Electron const* e) {
      return e->pt() > 20. && fabs(e->eta()) < 2.1 && e->isMVAwp80Electron();
    });

    boost::range::sort(electrons, DescendingPt);
    boost::range::sort(tracks, DescendingTrackPt);
    boost::range::sort(tracksforiso, DescendingTrackPt);

    ac::Candidate z_cand;
    if (electrons.size() >= 2) {
      z_cand.setVector(electrons[0]->vector() + electrons[1]->vector());
      z_cand.setCharge(electrons[0]->charge() + electrons[1]->charge());
    }

    if (electrons.size()==2 && veto_electrons.size() == 2 && z_cand.charge() == 0 && veto_muons.size()==0) {
      electron_overlaps_=0;
      muon_overlaps_=0;
      for(unsigned j=0; j<tracks.size(); j++){
        bool no_muon_overlap=1;
        for (unsigned i = 0; i<all_muons.size(); i++){
          no_muon_overlap&&(DeltaRTrack(tracks.at(j),all_muons.at(i))>0.3);
        }
        if (!no_muon_overlap) muon_overlaps_ +=1;
        bool no_electron_overlap=1;
        for (unsigned i = 0; i<all_electrons.size(); i++){
          no_electron_overlap&&(DeltaRTrack(tracks.at(j),all_electrons.at(i))>0.3);
        }
        if (!no_electron_overlap) electron_overlaps_ +=1;
      }

      pt_1_ = electrons[0]->pt();
      pt_2_ = electrons[1]->pt();
      eta_1_ = electrons[0]->eta();
      eta_2_ = electrons[1]->eta();
      m_ll_ = z_cand.M();
      pt_ll_ = z_cand.pt();
      dr_ll_ = DeltaR(electrons[0], electrons[1]);


      unsigned trg_lookup = is_data_ ? info->run() : year_;
      if (year_ == 2016) {
        auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele27_WPTight_Gsf");
        trg_1_ =
            IsFilterMatchedDR(electrons[0], trg_objs, filters_Ele27_.Lookup(trg_lookup), 0.3)&&electrons[0]->pt()>30;
        trg_2_ =
            IsFilterMatchedDR(electrons[1], trg_objs, filters_Ele27_.Lookup(trg_lookup), 0.3)&&electrons[1]->pt()>30;
      } else if (year_ == 2017) {
        auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele35_WPTight_Gsf");
        trg_1_ = IsFilterMatchedDR(electrons[0], trg_objs, filters_Ele35_.Lookup(trg_lookup), 0.3)&&electrons[0]->pt()>38;
        trg_2_ = IsFilterMatchedDR(electrons[1], trg_objs, filters_Ele35_.Lookup(trg_lookup), 0.3)&&electrons[1]->pt()>38;
      } else {
        auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_Ele32_WPTight_Gsf");
        trg_1_ = IsFilterMatchedDR(electrons[0], trg_objs, filters_Ele32_.Lookup(trg_lookup), 0.3)&&electrons[0]->pt()>35;
        trg_2_ = IsFilterMatchedDR(electrons[1], trg_objs, filters_Ele32_.Lookup(trg_lookup), 0.3)&&electrons[0]->pt()>35;
      }

      wt_pu_ = 1.;
      wt_1_ = 1.;
      wt_2_ = 1.;
      wt_trg1_ = 1.;
      wt_trg2_ = 1.;

      if (!is_data_) {
        auto const& pu_info = event->GetPtrVec<PileupInfo>("pileupInfo");
        for (PileupInfo const* pu : pu_info) {
          if (pu->bunchCrossing() == 0) {
            wt_pu_ = RooFunc(fns_["pileup_ratio"], {pu->trueNumInteractions()});
            break;
          }
        }
        wt_1_ = RooFunc(fns_["e_gsfidiso_ratio"], {pt_1_, eta_1_});
        wt_2_ = RooFunc(fns_["e_gsfidiso_ratio"], {pt_2_, eta_2_});
      }

      highestpt_pair_id_1_ = 0;
      highestpt_pair_id_2_ = 0;
      highestpt_pair_eta_ = -99;
      highestpt_pair_phi_ = -99;
      highestpt_pair_1_eta_ = -99;
      highestpt_pair_1_phi_ = -99;
      highestpt_pair_1_pt_ = -99;
      highestpt_pair_2_pt_ = -99;
      highestpt_pair_2_eta_ = -99;
      highestpt_pair_2_phi_ = -99;
      highestpt_pair_dR_=-99;
      highestpt_pair_mass_=-99;
      highestpt_pair_pt_=-99;
      highestpt_pair_reco_higgs_mass_=-99;
      highestpt_pair_reco_higgs_pt_=-99;


      std::vector<std::pair<std::pair<unsigned,unsigned>,double>> track_drs;
      std::vector<std::pair<std::pair<unsigned,unsigned>,double>> track_drs_smallcone;
      for (unsigned i = 0; i<tracks.size(); i++){
        if(DeltaRTrack(tracks.at(i),electrons.at(0))> 0.3 && DeltaRTrack(tracks.at(i),electrons.at(1))>0.3){
          for (unsigned j = 0; j < i; j++){
            if(DeltaRTrack(tracks.at(j),electrons.at(0))> 0.3 && DeltaRTrack(tracks.at(j),electrons.at(1))>0.3 && tracks.at(i)->charge()!=tracks.at(j)->charge() && (tracks.at(i)->pt()>10 || tracks.at(j)->pt()>10)){
              track_drs.push_back(std::make_pair(std::make_pair(i,j),DeltaRDiTrack(tracks.at(i),tracks.at(j))));
              if(DeltaRDiTrack(tracks.at(i),tracks.at(j))<0.1) track_drs_smallcone.push_back(std::make_pair(std::make_pair(i,j),(tracks.at(i)->vector()+tracks.at(j)->vector()).pt()));
            }
          } 
        }
      }
      boost::range::sort(track_drs, AscendingDR);
      boost::range::sort(track_drs_smallcone,DescendingPairPt);
      highestpt_pair_iso_ = -99; 
      highestpt_pair_looser_iso_ = -99; 
      wt_rhoiso_ = 1.;
      if(track_drs_smallcone.size()>0){
       highestpt_pair_id_1_=track_drs_smallcone.at(0).first.first;
       highestpt_pair_id_2_=track_drs_smallcone.at(0).first.second;
       highestpt_pair_reco_higgs_mass_=(electrons.at(0)->vector()+electrons.at(1)->vector()+tracks.at(highestpt_pair_id_1_)->vector()+tracks.at(highestpt_pair_id_2_)->vector()).M();
       highestpt_pair_reco_higgs_pt_=(electrons.at(0)->vector()+electrons.at(1)->vector()+tracks.at(highestpt_pair_id_1_)->vector()+tracks.at(highestpt_pair_id_2_)->vector()).pt();
       highestpt_pair_dR_=DeltaRDiTrack(tracks.at(track_drs_smallcone.at(0).first.first),tracks.at(track_drs_smallcone.at(0).first.second));
       highestpt_pair_mass_=(tracks.at(track_drs_smallcone.at(0).first.first)->vector()+tracks.at(track_drs_smallcone.at(0).first.second)->vector()).M();
       highestpt_pair_pt_=(tracks.at(track_drs_smallcone.at(0).first.first)->vector()+tracks.at(track_drs_smallcone.at(0).first.second)->vector()).pt();
       highestpt_pair_eta_=(tracks.at(track_drs_smallcone.at(0).first.first)->vector()+tracks.at(track_drs_smallcone.at(0).first.second)->vector()).eta();
       highestpt_pair_phi_=(tracks.at(track_drs_smallcone.at(0).first.first)->vector()+tracks.at(track_drs_smallcone.at(0).first.second)->vector()).phi();
       highestpt_pair_1_eta_=tracks.at(track_drs_smallcone.at(0).first.first)->eta();
       highestpt_pair_1_phi_=tracks.at(track_drs_smallcone.at(0).first.first)->phi();
       highestpt_pair_1_pt_=tracks.at(track_drs_smallcone.at(0).first.first)->pt();
       highestpt_pair_2_eta_=tracks.at(track_drs_smallcone.at(0).first.second)->eta();
       highestpt_pair_2_phi_=tracks.at(track_drs_smallcone.at(0).first.second)->phi();
       highestpt_pair_2_pt_=tracks.at(track_drs_smallcone.at(0).first.second)->pt();
       highestpt_pair_iso_=0;
       highestpt_pair_looser_iso_=0;
       Zrho_dphi_ = ROOT::Math::VectorUtil::DeltaPhi(tracks.at(track_drs_smallcone.at(0).first.first)->vector()+tracks.at(track_drs_smallcone.at(0).first.second)->vector(),z_cand.vector());
       for (unsigned i=0; i <tracks.size(); i++){
         if( i!=highestpt_pair_id_1_ && i!=highestpt_pair_id_2_){
           if(DeltaRTrackPair(tracks.at(highestpt_pair_id_1_),tracks.at(highestpt_pair_id_2_),tracks.at(i))<0.3) highestpt_pair_iso_+=tracks.at(i)->pt();
         }
       }
       for (unsigned i=0; i <tracksforiso.size(); i++){
         //if( i!=highestpt_pair_id_1_ && i!=highestpt_pair_id_2_){
           if(DeltaRTrackPair(tracks.at(highestpt_pair_id_1_),tracks.at(highestpt_pair_id_2_),tracksforiso.at(i))<0.3) highestpt_pair_looser_iso_+=tracksforiso.at(i)->pt();
         //}
       }
       highestpt_pair_looser_iso_-=highestpt_pair_1_pt_;
       highestpt_pair_looser_iso_-=highestpt_pair_2_pt_;
      }
      if(highestpt_pair_pt_>20.){
          wt_rhoiso_ = RooFunc(fns_["rhoiso_ratio_etainc"], {highestpt_pair_pt_, std::abs(highestpt_pair_eta_)});
      } else wt_rhoiso_ = RooFunc(fns_["rhoiso_ratio_etainc"], {21., std::abs(highestpt_pair_eta_)});

      tree_->Fill();
    }


    return 0;
  }
  int DiElectronMesonAnalysis::PostAnalysis() {
    return 0;
  }

  void DiElectronMesonAnalysis::PrintInfo() {}



}
