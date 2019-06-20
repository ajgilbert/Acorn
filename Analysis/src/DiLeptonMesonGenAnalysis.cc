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
#include "Acorn/Analysis/interface/DiLeptonMesonGenAnalysis.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/PileupInfo.h"
#include "Acorn/NTupler/interface/EventInfo.h"

namespace ac {

  DiLeptonMesonGenAnalysis::DiLeptonMesonGenAnalysis(std::string const& name)
      : ModuleBase(name), fs_(nullptr), year_(2016) {}

  DiLeptonMesonGenAnalysis::~DiLeptonMesonGenAnalysis() { ; }

  int DiLeptonMesonGenAnalysis::PreAnalysis() {
    if (fs_) {
      tree_ = fs_->make<TTree>("DiLeptonMesonGenAnalysis", "DiLeptonMesonGenAnalysis");
      tree_->Branch("pt_1", &pt_1_);
      tree_->Branch("pt_2", &pt_2_);
      tree_->Branch("eta_1", &eta_1_);
      tree_->Branch("eta_2", &eta_2_);
      tree_->Branch("m_ll", &m_ll_);
      tree_->Branch("pt_ll", &pt_ll_);
      tree_->Branch("dr_ll", &dr_ll_);
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
      tree_->Branch("reco_trk_looser_iso", &reco_trk_looser_iso_);
      tree_->Branch("reco_higgs_mass", &reco_higgs_mass_);
      tree_->Branch("reco_higgs_pt", &reco_higgs_pt_);
      tree_->Print();
    }

  return 0;
  }

  int DiLeptonMesonGenAnalysis::Execute(TreeEvent* event) {

    auto const* info = event->GetPtr<EventInfo>("eventInfo");

    std::vector<ac::Muon *> muons = event->GetPtrVec<ac::Muon>("muons");
    std::vector<ac::Electron *> electrons = event->GetPtrVec<ac::Electron>("electrons");
    std::vector<ac::Track *> tracks = event->GetPtrVec<ac::Track>("Tracks");
    std::vector<ac::Track *> tracksforiso = event->GetPtrVec<ac::Track>("TracksForIso");
    std::vector<GenParticle *> gen_parts = event->GetPtrVec<ac::GenParticle>("genParticles");
    std::vector<int> higgs_daughters;
    std::vector<int> rho_daughters;
    std::vector<GenParticle *> pions_from_meson;

    ac::keep_if(muons, [](ac::Muon const* m) {
      return m->pt() > 30. && fabs(m->eta()) < 2.4 && m->isMediumMuon() && MuonPFIso(m) < 0.15;
    });

    ac::keep_if(electrons, [](ac::Electron const* e) {
      return e->pt() > 30. && fabs(e->eta()) < 2.1 && e->isMVAwp80Electron();
    });

    boost::range::sort(muons, DescendingPt);
    boost::range::sort(electrons, DescendingPt);
    boost::range::sort(tracks, DescendingTrackPt);
    boost::range::sort(tracksforiso, DescendingTrackPt);

    ac::Candidate z_cand_mu;
    ac::Candidate z_cand_ele;
    if (muons.size() >= 2) {
      z_cand_mu.setVector(muons[0]->vector() + muons[1]->vector());
      z_cand_mu.setCharge(muons[0]->charge() + muons[1]->charge());
    }

    if (electrons.size() >= 2) {
      z_cand_ele.setVector(electrons[0]->vector() + electrons[1]->vector());
      z_cand_ele.setCharge(electrons[0]->charge() + electrons[1]->charge());
    }

    if ((muons.size() == 2 && z_cand_mu.charge() == 0 && electrons.size()==0) || (electrons.size() == 2 && z_cand_ele.charge() == 0 && muons.size()==0)) {
      if (muons.size()==2){
          pt_1_ = muons[0]->pt();
          pt_2_ = muons[1]->pt();
          eta_1_ = muons[0]->eta();
          eta_2_ = muons[1]->eta();
          m_ll_ = z_cand_mu.M();
          pt_ll_ = z_cand_mu.pt();
          dr_ll_ = DeltaR(muons[0], muons[1]);
       } else {
          pt_1_ = electrons[0]->pt();
          pt_2_ = electrons[1]->pt();
          eta_1_ = electrons[0]->eta();
          eta_2_ = electrons[1]->eta();
          m_ll_ = z_cand_ele.M();
          pt_ll_ = z_cand_ele.pt();
          dr_ll_ = DeltaR(electrons[0], electrons[1]);
       }


       for (auto const& part : gen_parts) {
          if(part->pdgId()==25&&part->statusFlags().isLastCopy()){
              higgs_daughters = part->daughters();
              for ( unsigned int i = 0; i < higgs_daughters.size(); i++){
                if( abs(gen_parts.at(higgs_daughters.at(i))->pdgId())== 113 || abs(gen_parts.at(higgs_daughters.at(i))->pdgId())==223 || abs(gen_parts.at(higgs_daughters.at(i))->pdgId())==333){
                  rho_daughters = gen_parts.at(higgs_daughters.at(i))->daughters();
                  for ( unsigned int j =0; j < rho_daughters.size(); j++){
                    if(abs(gen_parts.at(rho_daughters.at(j))->pdgId())==211 || abs(gen_parts.at(rho_daughters.at(j))->pdgId())==321) pions_from_meson.push_back(gen_parts.at(rho_daughters.at(j)));
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
      reco_higgs_mass_=-99;
      reco_higgs_pt_=-99;


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
            auto v1 = ROOT::Math::PtEtaPhiMVector(tracks.at(i)->vector());
            v1.SetM(0.493);
            auto v2 = ROOT::Math::PtEtaPhiMVector(tracks.at(reco_trk_id_1_)->vector());
            v2.SetM(0.493);
            reco_trk_mass_=(v1+v2).M();
            //reco_trk_mass_=((tracks.at(i)->vector()).SetM(0.493)+(tracks.at(reco_trk_id_1_)->vector()).SetM(0.493)).M();
            reco_trk_pt_=(tracks.at(i)->vector()+tracks.at(reco_trk_id_1_)->vector()).pt();
            if(muons.size()==2){
              reco_higgs_mass_=(muons.at(0)->vector()+muons.at(1)->vector()+tracks.at(i)->vector()+tracks.at(reco_trk_id_1_)->vector()).M();
              reco_higgs_pt_=(muons.at(0)->vector()+muons.at(1)->vector()+tracks.at(i)->vector()+tracks.at(reco_trk_id_1_)->vector()).pt();
            } else {
              reco_higgs_mass_=(electrons.at(0)->vector()+electrons.at(1)->vector()+tracks.at(i)->vector()+tracks.at(reco_trk_id_1_)->vector()).M();
              reco_higgs_pt_=(electrons.at(0)->vector()+electrons.at(1)->vector()+tracks.at(i)->vector()+tracks.at(reco_trk_id_1_)->vector()).pt();
           }
         }
       }
      }
      reco_trk_iso_=-99;
      reco_trk_looser_iso_=-99;
      if(reco_trk_pt_1_>-1 && reco_trk_pt_2_>-1){
        reco_trk_iso_=0;
        reco_trk_looser_iso_=0;
        for ( unsigned i=0; i<tracks.size() ; i++){
          if(i!=reco_trk_id_1_ && i!=reco_trk_id_2_){
            if(DeltaRTrackPair(tracks.at(reco_trk_id_1_),tracks.at(reco_trk_id_2_),tracks.at(i))<0.3) reco_trk_iso_+=tracks.at(i)->pt();
           }
        }
        for ( unsigned i=0; i<tracksforiso.size() ; i++){
            if(DeltaRTrackPair(tracks.at(reco_trk_id_1_),tracks.at(reco_trk_id_2_),tracksforiso.at(i))<0.3) reco_trk_looser_iso_+=tracksforiso.at(i)->pt();
        }
       reco_trk_looser_iso_-=reco_trk_pt_1_;
       reco_trk_looser_iso_-=reco_trk_pt_2_;
      }
      tree_->Fill();
    }


    return 0;
  }
  int DiLeptonMesonGenAnalysis::PostAnalysis() {
    return 0;
  }

  void DiLeptonMesonGenAnalysis::PrintInfo() {}



}
