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
#include "Acorn/Analysis/interface/DiMuonMesonAnalysis.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"

namespace ac {

DiMuonMesonAnalysis::DiMuonMesonAnalysis(std::string const& name)
    : ModuleBase(name), fs_(nullptr), year_(2016), is_data_(true) {}

DiMuonMesonAnalysis::~DiMuonMesonAnalysis() { ; }

int DiMuonMesonAnalysis::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("DiMuonMesonAnalysis", "DiMuonMesonAnalysis");
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
    tree_->Branch("gen_trk_pt_1", &gen_trk_pt_1_);
    tree_->Branch("gen_trk_pt_2", &gen_trk_pt_2_);
    tree_->Branch("gen_trk_eta_1", &gen_trk_eta_1_);
    tree_->Branch("gen_trk_eta_2", &gen_trk_eta_2_);
    tree_->Branch("gen_trk_phi_1", &gen_trk_phi_1_);
    tree_->Branch("gen_trk_phi_2", &gen_trk_phi_2_);
    tree_->Branch("reco_trk_pt_1", &reco_trk_pt_1_);
    tree_->Branch("reco_trk_pt_2", &reco_trk_pt_2_);
    tree_->Branch("reco_trk_eta_1", &reco_trk_eta_1_);
    tree_->Branch("reco_trk_eta_2", &reco_trk_eta_2_);
    tree_->Branch("reco_trk_phi_1", &reco_trk_phi_1_);
    tree_->Branch("reco_trk_phi_2", &reco_trk_phi_2_);
    tree_->Branch("reco_trk_charge_1", &reco_trk_charge_1_);
    tree_->Branch("reco_trk_charge_2", &reco_trk_charge_2_);
    tree_->Branch("reco_trk_id_1", &reco_trk_id_1_);
    tree_->Branch("reco_trk_id_2", &reco_trk_id_2_);
    tree_->Branch("reco_trk_dR", &reco_trk_dR_);
    tree_->Branch("reco_trk_mass", &reco_trk_mass_);
    tree_->Branch("reco_trk_pt", &reco_trk_pt_);
    tree_->Branch("reco_trk_iso", &reco_trk_iso_);
    tree_->Branch("closest_pair_id_1", &closest_pair_id_1_);
    tree_->Branch("closest_pair_id_2", &closest_pair_id_2_);
    tree_->Branch("closest_pair_dR", &closest_pair_dR_);
    tree_->Branch("closest_pair_mass", &closest_pair_mass_);
    tree_->Branch("closest_pair_pt", &closest_pair_pt_);
    tree_->Branch("nclosest_pair_id_1", &nclosest_pair_id_1_);
    tree_->Branch("nclosest_pair_id_2", &nclosest_pair_id_2_);
    tree_->Branch("nclosest_pair_dR", &nclosest_pair_dR_);
    tree_->Branch("nclosest_pair_mass", &nclosest_pair_mass_);
    tree_->Branch("nclosest_pair_pt", &nclosest_pair_pt_);
    tree_->Branch("nnclosest_pair_id_1", &nnclosest_pair_id_1_);
    tree_->Branch("nnclosest_pair_id_2", &nnclosest_pair_id_2_);
    tree_->Branch("nnclosest_pair_dR", &nnclosest_pair_dR_);
    tree_->Branch("nnclosest_pair_mass", &nnclosest_pair_mass_);
    tree_->Branch("nnclosest_pair_pt", &nnclosest_pair_pt_);
    tree_->Branch("highestpt_pair_id_1", &highestpt_pair_id_1_);
    tree_->Branch("highestpt_pair_id_2", &highestpt_pair_id_2_);
    tree_->Branch("highestpt_pair_dR", &highestpt_pair_dR_);
    tree_->Branch("highestpt_pair_mass", &highestpt_pair_mass_);
    tree_->Branch("highestpt_pair_pt", &highestpt_pair_pt_);
    tree_->Branch("highestpt_pair_iso", &highestpt_pair_iso_);
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

  int DiMuonMesonAnalysis::Execute(TreeEvent* event) {

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    std::vector<ac::Muon *> muons = event->GetPtrVec<ac::Muon>("muons");
    std::vector<ac::Track *> tracks = event->GetPtrVec<ac::Track>("Tracks");
    std::vector<GenParticle *> gen_parts = event->GetPtrVec<ac::GenParticle>("genParticles");
    std::vector<int> higgs_daughters;
    std::vector<int> rho_daughters;
    std::vector<GenParticle *> pions_from_meson;

    // Apply pT and ID cuts
    // The medium ID already includes iso?
    ac::keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && MuonPFIso(m) < 0.15;
    });

    boost::range::sort(muons, DescendingPt);
    boost::range::sort(tracks, DescendingTrackPt);

    ac::Candidate z_cand;
    if (muons.size() >= 2) {
      z_cand.setVector(muons[0]->vector() + muons[1]->vector());
      z_cand.setCharge(muons[0]->charge() + muons[1]->charge());
    }

    if (muons.size() == 2 && z_cand.charge() == 0) {
      pt_1_ = muons[0]->pt();
      pt_2_ = muons[1]->pt();
      eta_1_ = muons[0]->eta();
      eta_2_ = muons[1]->eta();
      m_ll_ = z_cand.M();
      pt_ll_ = z_cand.pt();
      dr_ll_ = DeltaR(muons[0], muons[1]);


      unsigned trg_lookup = is_data_ ? info->run() : year_;
      if (year_ == 2016) {
        auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu24");
        auto const& trg_objs_tk = event->GetPtrVec<TriggerObject>("triggerObjects_IsoTkMu24");
        trg_1_ =
            IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3) ||
            IsFilterMatchedDR(muons[0], trg_objs_tk, filters_IsoTkMu24_.Lookup(trg_lookup), 0.3);
        trg_2_ =
            IsFilterMatchedDR(muons[1], trg_objs, filters_IsoMu24_.Lookup(trg_lookup), 0.3) ||
            IsFilterMatchedDR(muons[1], trg_objs_tk, filters_IsoTkMu24_.Lookup(trg_lookup), 0.3);
      } else {
        auto const& trg_objs = event->GetPtrVec<TriggerObject>("triggerObjects_IsoMu27");
        trg_1_ = IsFilterMatchedDR(muons[0], trg_objs, filters_IsoMu27_.Lookup(trg_lookup), 0.3);
        trg_2_ = IsFilterMatchedDR(muons[1], trg_objs, filters_IsoMu27_.Lookup(trg_lookup), 0.3);
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
        wt_1_ = RooFunc(fns_["m_idisotrk_ratio"], {pt_1_, eta_1_});
        wt_2_ = RooFunc(fns_["m_idisotrk_ratio"], {pt_2_, eta_2_});
        wt_trg1_ = RooFunc(fns_["m_trg_ratio"], {pt_1_, eta_1_});
        wt_trg2_ = RooFunc(fns_["m_trg_ratio"], {pt_2_, eta_2_});

       for (auto const& part : gen_parts) {
          if(part->pdgId()==25&&part->statusFlags().isLastCopy()){
              higgs_daughters = part->daughters();
              for ( unsigned int i = 0; i < higgs_daughters.size(); i++){
                if( abs(gen_parts.at(higgs_daughters.at(i))->pdgId())== 113 || abs(gen_parts.at(higgs_daughters.at(i))->pdgId())==223){
                  rho_daughters = gen_parts.at(higgs_daughters.at(i))->daughters();
                  for ( unsigned int j =0; j < rho_daughters.size(); j++){
                    if(abs(gen_parts.at(rho_daughters.at(j))->pdgId())==211) pions_from_meson.push_back(gen_parts.at(rho_daughters.at(j)));
                  }
               }
            }
          }
        }
       boost::range::sort(pions_from_meson, DescendingPt);
       if(pions_from_meson.size() > 1){
          gen_trk_pt_1_=pions_from_meson.at(0)->pt();
          gen_trk_eta_1_=pions_from_meson.at(0)->eta();
          gen_trk_phi_1_=pions_from_meson.at(0)->phi();
          gen_trk_pt_2_=pions_from_meson.at(1)->pt();
          gen_trk_eta_2_=pions_from_meson.at(1)->eta();
          gen_trk_phi_2_=pions_from_meson.at(1)->phi();
        } else {
          gen_trk_pt_1_=-1;
          gen_trk_eta_1_=-1;
          gen_trk_phi_1_=-1;
          gen_trk_pt_2_=-1;
          gen_trk_eta_2_=-1;
          gen_trk_phi_2_=-1;
       }
      }

      reco_trk_pt_1_=-99;
      reco_trk_eta_1_=-99;
      reco_trk_phi_1_=-99;
      reco_trk_charge_1_=-99;
      reco_trk_id_1_=0;
      reco_trk_pt_2_=-99;
      reco_trk_eta_2_=-99;
      reco_trk_phi_2_=-99;
      reco_trk_charge_2_=-99;
      reco_trk_id_2_=0;
      reco_trk_dR_=-99;
      reco_trk_mass_=-99;
      reco_trk_pt_=-99;
      highestpt_pair_id_1_ = 0;
      highestpt_pair_id_2_ = 0;
      highestpt_pair_dR_=-99;
      highestpt_pair_mass_=-99;
      highestpt_pair_pt_=-99;

      for (unsigned i = 0; i<tracks.size(); i++){
        double mindr=100.;
        double thedr=0.;
        for (unsigned j = 0; j< pions_from_meson.size(); j++){
          thedr = DeltaRTrack(tracks.at(i),pions_from_meson.at(j));
          if( thedr < mindr) mindr=thedr;
        }
        if(mindr<0.01){
          if(reco_trk_pt_1_==-99){
            reco_trk_pt_1_ = tracks.at(i)->pt();
            reco_trk_eta_1_ = tracks.at(i)->eta();
            reco_trk_phi_1_ = tracks.at(i)->phi();
            reco_trk_charge_1_ = tracks.at(i)->charge();
            reco_trk_id_1_ = i;
          } else {
            reco_trk_pt_2_ = tracks.at(i)->pt();
            reco_trk_eta_2_ = tracks.at(i)->eta();
            reco_trk_phi_2_ = tracks.at(i)->phi();
            reco_trk_charge_2_ = tracks.at(i)->charge();
            reco_trk_id_2_ = i;
            reco_trk_dR_=DeltaRDiTrack(tracks.at(i),tracks.at(reco_trk_id_1_));
            reco_trk_mass_=(tracks.at(i)->vector()+tracks.at(reco_trk_id_1_)->vector()).M();
            reco_trk_pt_=(tracks.at(i)->vector()+tracks.at(reco_trk_id_1_)->vector()).pt();
         }
       }
      }
      reco_trk_iso_=-99;
      if(reco_trk_pt_1_>-1 && reco_trk_pt_2_>-1){
        reco_trk_iso_=0;
        for ( unsigned i=0; i<tracks.size() ; i++){
          if(i!=reco_trk_id_1_ && i!=reco_trk_id_2_){
            if(DeltaRTrackPair(tracks.at(reco_trk_id_1_),tracks.at(reco_trk_id_2_),tracks.at(i))<0.3) reco_trk_iso_+=tracks.at(i)->pt();
           }
        }
      }
      std::vector<std::pair<std::pair<unsigned,unsigned>,double>> track_drs;
      std::vector<std::pair<std::pair<unsigned,unsigned>,double>> track_drs_smallcone;
      for (unsigned i = 0; i<tracks.size(); i++){
        if(DeltaRTrack(tracks.at(i),muons.at(0))> 0.3 && DeltaRTrack(tracks.at(i),muons.at(1))>0.3){
          for (unsigned j = 0; j < i; j++){
            if(DeltaRTrack(tracks.at(j),muons.at(0))> 0.3 && DeltaRTrack(tracks.at(j),muons.at(1))>0.3 && tracks.at(i)->charge()!=tracks.at(j)->charge()){
              track_drs.push_back(std::make_pair(std::make_pair(i,j),DeltaRDiTrack(tracks.at(i),tracks.at(j))));
              if(DeltaRDiTrack(tracks.at(i),tracks.at(j))<0.1) track_drs_smallcone.push_back(std::make_pair(std::make_pair(i,j),(tracks.at(i)->vector()+tracks.at(j)->vector()).pt()));
            }
          } 
        }
      }
      boost::range::sort(track_drs, AscendingDR);
      boost::range::sort(track_drs_smallcone,DescendingPairPt);
      highestpt_pair_iso_ = -99; 
      if(track_drs_smallcone.size()>0){
       highestpt_pair_id_1_=track_drs_smallcone.at(0).first.first;
       highestpt_pair_id_2_=track_drs_smallcone.at(0).first.second;
       highestpt_pair_dR_=DeltaRDiTrack(tracks.at(track_drs_smallcone.at(0).first.first),tracks.at(track_drs_smallcone.at(0).first.second));
       highestpt_pair_mass_=(tracks.at(track_drs_smallcone.at(0).first.first)->vector()+tracks.at(track_drs_smallcone.at(0).first.second)->vector()).M();
       highestpt_pair_pt_=(tracks.at(track_drs_smallcone.at(0).first.first)->vector()+tracks.at(track_drs_smallcone.at(0).first.second)->vector()).pt();
       highestpt_pair_iso_=0;
       for (unsigned i=0; i <tracks.size(); i++){
         if( i!=highestpt_pair_id_1_ && i!=highestpt_pair_id_2_){
           if(DeltaRTrackPair(tracks.at(highestpt_pair_id_1_),tracks.at(highestpt_pair_id_2_),tracks.at(i))<0.3) highestpt_pair_iso_+=tracks.at(i)->pt();
         }
       }
      }
      if(track_drs.size()>2){
        closest_pair_id_1_=track_drs.at(0).first.first;
        closest_pair_id_2_=track_drs.at(0).first.second;
        closest_pair_dR_ = track_drs.at(0).second;
        closest_pair_mass_ = (tracks.at(track_drs.at(0).first.first)->vector()+tracks.at(track_drs.at(0).first.second)->vector()).M();
        closest_pair_pt_ = (tracks.at(track_drs.at(0).first.first)->vector()+tracks.at(track_drs.at(0).first.second)->vector()).pt();
        nclosest_pair_id_1_=track_drs.at(1).first.first;
        nclosest_pair_id_2_=track_drs.at(1).first.second;
        nclosest_pair_dR_ = track_drs.at(1).second;
        nclosest_pair_mass_ = (tracks.at(track_drs.at(1).first.first)->vector()+tracks.at(track_drs.at(1).first.second)->vector()).M();
        nclosest_pair_pt_ = (tracks.at(track_drs.at(1).first.first)->vector()+tracks.at(track_drs.at(1).first.second)->vector()).pt();
        nnclosest_pair_id_1_=track_drs.at(2).first.first;
        nnclosest_pair_id_2_=track_drs.at(2).first.second;
        nnclosest_pair_dR_ = track_drs.at(2).second;
        nnclosest_pair_mass_ = (tracks.at(track_drs.at(2).first.first)->vector()+tracks.at(track_drs.at(2).first.second)->vector()).M();
        nnclosest_pair_pt_ = (tracks.at(track_drs.at(2).first.first)->vector()+tracks.at(track_drs.at(2).first.second)->vector()).pt();
      } else {
        closest_pair_id_1_=0;
        closest_pair_id_2_=0;
        closest_pair_dR_ = 0;
        closest_pair_mass_ = 0;
        closest_pair_pt_ = 0;
        nclosest_pair_id_1_=0;
        nclosest_pair_id_2_=0;
        nclosest_pair_dR_ = 0;
        nclosest_pair_mass_ = 0;
        nclosest_pair_pt_ = 0;
        nnclosest_pair_id_1_=0;
        nnclosest_pair_id_2_=0;
        nnclosest_pair_dR_ = 0;
        nnclosest_pair_mass_ = 0;
        nnclosest_pair_pt_ = 0;
      }
      tree_->Fill();
    }


    return 0;
  }
  int DiMuonMesonAnalysis::PostAnalysis() {
    return 0;
  }

  void DiMuonMesonAnalysis::PrintInfo() {}



}
