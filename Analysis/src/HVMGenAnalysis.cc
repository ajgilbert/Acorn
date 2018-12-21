#include <algorithm>
#include <map>
#include "TMath.h"
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"
#include "RooRealVar.h"
#include "Acorn/Analysis/interface/HVMGenAnalysis.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

namespace ac {

  HVMGenAnalysis::HVMGenAnalysis(std::string const& name)
      : ModuleBase(name), fs_(nullptr){}

  HVMGenAnalysis::~HVMGenAnalysis() { ; }

  int HVMGenAnalysis::PreAnalysis() {
    if (fs_) {
      tree_ = fs_->make<TTree>("HVMGenAnalysis", "HVMGenAnalysis");
      tree_->Branch("muon1_pt", &muon1_pt_);
      tree_->Branch("muon2_pt", &muon2_pt_);
      tree_->Branch("muon1_eta", &muon1_eta_);
      tree_->Branch("muon2_eta", &muon2_eta_);
      tree_->Branch("muon1_phi", &muon1_phi_);
      tree_->Branch("muon2_phi", &muon2_phi_);
      tree_->Branch("pion1_pt", &pion1_pt_);
      tree_->Branch("pion2_pt", &pion2_pt_);
      tree_->Branch("pion3_pt", &pion3_pt_);
      tree_->Branch("pion4_pt", &pion4_pt_);
      tree_->Branch("pion5_pt", &pion5_pt_);
      tree_->Branch("pion1_eta", &pion1_eta_);
      tree_->Branch("pion2_eta", &pion2_eta_);
      tree_->Branch("pion3_eta", &pion3_eta_);
      tree_->Branch("pion4_eta", &pion4_eta_);
      tree_->Branch("pion5_eta", &pion5_eta_);
      tree_->Branch("pion1_phi", &pion1_phi_);
      tree_->Branch("pion2_phi", &pion2_phi_);
      tree_->Branch("pion3_phi", &pion3_phi_);
      tree_->Branch("pion4_phi", &pion4_phi_);
      tree_->Branch("pion5_phi", &pion5_phi_);
      tree_->Branch("pion1_charge", &pion1_charge_);
      tree_->Branch("pion2_charge", &pion2_charge_);
      tree_->Branch("pion3_charge", &pion3_charge_);
      tree_->Branch("pion4_charge", &pion4_charge_);
      tree_->Branch("pion5_charge", &pion5_charge_);
      tree_->Branch("pion1_fromres", &pion1_fromres_);
      tree_->Branch("pion2_fromres", &pion2_fromres_);
      tree_->Branch("pion3_fromres", &pion3_fromres_);
      tree_->Branch("pion4_fromres", &pion4_fromres_);
      tree_->Branch("pion5_fromres", &pion5_fromres_);
      tree_->Branch("res_mass", &res_mass_);
      tree_->Branch("res_pt", &res_pt_);
      tree_->Branch("res_eta", &res_eta_);
      tree_->Branch("res_phi", &res_phi_);
      tree_->Branch("genres_mass", &genres_mass_);
      tree_->Branch("genres_pt", &genres_pt_);
      tree_->Branch("genres_eta", &genres_eta_);
      tree_->Branch("genres_phi", &genres_phi_);
      tree_->Branch("pions_dr", &pions_dr_);

      // tree_->Branch("wt", &wt_);
      // tree_->Branch("n_jets", &n_jets_);
    }
    return 0;
  }

  int HVMGenAnalysis::PostAnalysis() {
    return 0;
  }

  int HVMGenAnalysis::Execute(TreeEvent* event) {
    res_mass_=-999;
    res_pt_=-999;
    res_eta_=-999;
    res_phi_=-999;
    muon1_pt_=-999;
    muon2_pt_=-999;
    muon1_eta_=-999;
    muon2_eta_=-999;
    muon1_phi_=-999;
    muon2_phi_=-999;
    pion1_pt_=-999;
    pion2_pt_=-999;
    pion3_pt_=-999;
    pion4_pt_=-999;
    pion5_pt_=-999;
    pion1_eta_=-999;
    pion2_eta_=-999;
    pion3_eta_=-999;
    pion4_eta_=-999;
    pion5_eta_=-999;
    pion1_phi_=-999;
    pion2_phi_=-999;
    pion3_phi_=-999;
    pion4_phi_=-999;
    pion5_phi_=-999;
    pion1_charge_=-999;
    pion2_charge_=-999;
    pion3_charge_=-999;
    pion4_charge_=-999;
    pion5_charge_=-999;
    pion1_fromres_=0;
    pion2_fromres_=0;
    pion3_fromres_=0;
    pion4_fromres_=0;
    pion5_fromres_=0;


    std::vector<GenParticle *> gen_parts;
    std::vector<GenParticle *> muons;
    std::vector<GenParticle *> pions;
    std::vector<int> higgs_daughters;
    std::vector<int> rho_daughters;
    std::vector<GenParticle *> pions_from_meson;
    int meson_index = 0;
    gen_parts = event->GetPtrVec<ac::GenParticle>("genParticles");
    for (auto const& part : gen_parts) {
      if(part->pdgId()==25&&part->statusFlags().isLastCopy()){
          higgs_daughters = part->daughters();
          for ( unsigned int i = 0; i < higgs_daughters.size(); i++){
            if( abs(gen_parts.at(higgs_daughters.at(i))->pdgId())== 113 || abs(gen_parts.at(higgs_daughters.at(i))->pdgId())==223){
              meson_index = higgs_daughters.at(i);
              rho_daughters = gen_parts.at(higgs_daughters.at(i))->daughters();
              genres_mass_ = gen_parts.at(higgs_daughters.at(i))->M();
              genres_pt_ = gen_parts.at(higgs_daughters.at(i))->pt();
              genres_eta_ = gen_parts.at(higgs_daughters.at(i))->eta();
              genres_phi_ = gen_parts.at(higgs_daughters.at(i))->phi();
              for ( unsigned int j =0; j < rho_daughters.size(); j++){
                if(abs(gen_parts.at(rho_daughters.at(j))->pdgId())==211) pions_from_meson.push_back(gen_parts.at(rho_daughters.at(j)));
              }
           }
        }
      }
      if(part->status()==1){
        if(abs(part->pdgId())==13 && part->statusFlags().isPrompt()){
            muons.push_back(part);
        }
        if(abs(part->pdgId())==211 && part->statusFlags().isDirectHadronDecayProduct()){
            pions.push_back(part);
        }
      }
   }
   std::sort(muons.begin(), muons.end(), [&](ac::GenParticle const* c1, ac::GenParticle const* c2) {
      return c1->pt() > c2->pt();
    });
   std::sort(pions.begin(), pions.end(), [&](ac::GenParticle const* c1, ac::GenParticle const* c2) {
      return c1->pt() > c2->pt();
    });

    std::vector<int> mothers1_;
    std::vector<int> mothers2_;

  if(muons.size()>2){
    std::cout<<"There are "<<muons.size()<<" muons in this event! "<<std::endl;
    muon1_pt_=muons.at(0)->pt();
    muon2_pt_=muons.at(1)->pt();
    muon1_phi_=muons.at(0)->phi();
    muon2_phi_=muons.at(1)->phi();
    muon1_eta_=muons.at(0)->eta();
    muon2_eta_=muons.at(1)->eta();
  } else if (muons.size()==2){
    muon1_pt_=muons.at(0)->pt();
    muon2_pt_=muons.at(1)->pt();
    muon1_phi_=muons.at(0)->phi();
    muon2_phi_=muons.at(1)->phi();
    muon1_eta_=muons.at(0)->eta();
    muon2_eta_=muons.at(1)->eta();
 } else if (muons.size()==1){
    muon1_pt_=muons.at(0)->pt();
    muon1_phi_=muons.at(0)->phi();
    muon1_eta_=muons.at(0)->eta();
 }



 std::vector<int> mothers1;
 std::vector<int> mothers2;
 std::vector<int> mothers3;
 std::vector<int> mothers4;
 std::vector<int> mothers5;

 if(pions.size()>5){
    mothers1=pions.at(0)->mothers();
    for ( unsigned int i = 0 ; i < mothers1.size() ; i++){
      if( mothers1.at(i)==meson_index) pion1_fromres_ = 1;
    }
    mothers2=pions.at(1)->mothers();
    for ( unsigned int i = 0 ; i < mothers2.size() ; i++){
      if( mothers2.at(i)==meson_index) pion2_fromres_ = 1;
    }
    mothers3=pions.at(2)->mothers();
    for ( unsigned int i = 0 ; i < mothers3.size() ; i++){
      if( mothers3.at(i) ==meson_index) pion3_fromres_ = 1;
    }
    mothers4=pions.at(3)->mothers();
    for ( unsigned int i = 0 ; i < mothers4.size() ; i++){
      if(mothers4.at(i) == meson_index) pion4_fromres_ = 1;
    }
    mothers5=pions.at(4)->mothers();
    for ( unsigned int i = 0 ; i < mothers5.size() ; i++){
      if(mothers5.at(i) == meson_index) pion5_fromres_ = 1;
    }

    pion1_pt_=pions.at(0)->pt();
    pion1_eta_=pions.at(0)->eta();
    pion1_phi_=pions.at(0)->phi();
    pion1_charge_=pions.at(0)->charge();

    pion2_pt_=pions.at(1)->pt();
    pion2_eta_=pions.at(1)->eta();
    pion2_phi_=pions.at(1)->phi();
    pion2_charge_=pions.at(1)->charge();
    pion3_pt_=pions.at(2)->pt();
    pion3_eta_=pions.at(2)->eta();
    pion3_phi_=pions.at(2)->phi();
    pion3_charge_=pions.at(2)->charge();
    pion4_pt_=pions.at(3)->pt();
    pion4_eta_=pions.at(3)->eta();
    pion4_phi_=pions.at(3)->phi();
    pion4_charge_=pions.at(3)->charge();
    pion5_pt_=pions.at(4)->pt();
    pion5_eta_=pions.at(4)->eta();
    pion5_phi_=pions.at(4)->phi();
    pion5_charge_=pions.at(4)->charge();
 } else if (pions.size()==5){
    mothers1=pions.at(0)->mothers();
    for ( unsigned int i = 0 ; i < mothers1.size() ; i++){
      if(mothers1.at(i)==meson_index) pion1_fromres_ = 1;
    }
    mothers2=pions.at(1)->mothers();
    for ( unsigned int i = 0 ; i < mothers2.size() ; i++){
      if(mothers2.at(i)==meson_index) pion2_fromres_ = 1;
    }
    mothers3=pions.at(2)->mothers();
    for ( unsigned int i = 0 ; i < mothers3.size() ; i++){
      if(mothers3.at(i)==meson_index) pion3_fromres_ = 1;
    }
    mothers4=pions.at(3)->mothers();
    for ( unsigned int i = 0 ; i < mothers4.size() ; i++){
      if(mothers4.at(i)==meson_index) pion4_fromres_ = 1;
    }
    mothers5=pions.at(4)->mothers();
    for ( unsigned int i = 0 ; i < mothers5.size() ; i++){
      if(mothers5.at(i)==meson_index) pion5_fromres_ = 1;
    }

    pion1_pt_=pions.at(0)->pt();
    pion1_eta_=pions.at(0)->eta();
    pion1_phi_=pions.at(0)->phi();
    pion1_charge_=pions.at(0)->charge();
    pion2_pt_=pions.at(1)->pt();
    pion2_eta_=pions.at(1)->eta();
    pion2_phi_=pions.at(1)->phi();
    pion2_charge_=pions.at(1)->charge();
    pion3_pt_=pions.at(2)->pt();
    pion3_eta_=pions.at(2)->eta();
    pion3_phi_=pions.at(2)->phi();
    pion3_charge_=pions.at(2)->charge();
    pion4_pt_=pions.at(3)->pt();
    pion4_eta_=pions.at(3)->eta();
    pion4_phi_=pions.at(3)->phi();
    pion4_charge_=pions.at(3)->charge();
    pion5_pt_=pions.at(4)->pt();
    pion5_eta_=pions.at(4)->eta();
    pion5_phi_=pions.at(4)->phi();
    pion5_charge_=pions.at(4)->charge();
 } else if (pions.size()==4){
    mothers1=pions.at(0)->mothers();
    for ( unsigned int i = 0 ; i < mothers1.size() ; i++){
      if(mothers1.at(i)==meson_index) pion1_fromres_ = 1;
    }
    mothers2=pions.at(1)->mothers();
    for ( unsigned int i = 0 ; i < mothers2.size() ; i++){
      if(mothers2.at(i)==meson_index) pion2_fromres_ = 1;
    }
    mothers3=pions.at(2)->mothers();
    for ( unsigned int i = 0 ; i < mothers3.size() ; i++){
      if(mothers3.at(i)==meson_index) pion3_fromres_ = 1;
    }
    mothers4=pions.at(3)->mothers();
    for ( unsigned int i = 0 ; i < mothers4.size() ; i++){
      if(mothers4.at(i)==meson_index) pion4_fromres_ = 1;
    }

    pion1_pt_=pions.at(0)->pt();
    pion1_eta_=pions.at(0)->eta();
    pion1_phi_=pions.at(0)->phi();
    pion1_charge_=pions.at(0)->charge();
    pion2_pt_=pions.at(1)->pt();
    pion2_eta_=pions.at(1)->eta();
    pion2_phi_=pions.at(1)->phi();
    pion2_charge_=pions.at(1)->charge();
    pion3_pt_=pions.at(2)->pt();
    pion3_eta_=pions.at(2)->eta();
    pion3_phi_=pions.at(2)->phi();
    pion3_charge_=pions.at(2)->charge();
    pion4_pt_=pions.at(3)->pt();
    pion4_eta_=pions.at(3)->eta();
    pion4_phi_=pions.at(3)->phi();
    pion4_charge_=pions.at(3)->charge();
 } else if (pions.size()==3){
    mothers1=pions.at(0)->mothers();
    for ( unsigned int i = 0 ; i < mothers1.size() ; i++){
      if(mothers1.at(i)==meson_index) pion1_fromres_ = 1;
    }
    mothers2=pions.at(1)->mothers();
    for ( unsigned int i = 0 ; i < mothers2.size() ; i++){
      if(mothers2.at(i)==meson_index) pion2_fromres_ = 1;
    }
    mothers3=pions.at(2)->mothers();
    for ( unsigned int i = 0 ; i < mothers3.size() ; i++){
      if(mothers3.at(i)==meson_index) pion3_fromres_ = 1;
    }

    pion1_pt_=pions.at(0)->pt();
    pion1_eta_=pions.at(0)->eta();
    pion1_phi_=pions.at(0)->phi();
    pion1_charge_=pions.at(0)->charge();
    pion2_pt_=pions.at(1)->pt();
    pion2_eta_=pions.at(1)->eta();
    pion2_phi_=pions.at(1)->phi();
    pion2_charge_=pions.at(1)->charge();
    pion3_pt_=pions.at(2)->pt();
    pion3_eta_=pions.at(2)->eta();
    pion3_phi_=pions.at(2)->phi();
    pion3_charge_=pions.at(2)->charge();
 } else if (pions.size()==2){
    mothers1=pions.at(0)->mothers();
    for ( unsigned int i = 0 ; i < mothers1.size() ; i++){
      if(mothers1.at(i)==meson_index) pion1_fromres_ = 1;
    }
    mothers2=pions.at(1)->mothers();
    for ( unsigned int i = 0 ; i < mothers2.size() ; i++){
      if(mothers2.at(i)==meson_index) pion2_fromres_ = 1;
    }
    pion1_pt_=pions.at(0)->pt();
    pion1_eta_=pions.at(0)->eta();
    pion1_phi_=pions.at(0)->phi();
    pion1_charge_=pions.at(0)->charge();
    pion2_pt_=pions.at(1)->pt();
    pion2_eta_=pions.at(1)->eta();
    pion2_phi_=pions.at(1)->phi();
    pion2_charge_=pions.at(1)->charge();
  } else if (pions.size()==1){
    mothers1=pions.at(0)->mothers();
    for ( unsigned int i = 0 ; i < mothers1.size() ; i++){
      if(mothers1.at(i)==meson_index) pion1_fromres_ = 1;
    }
    pion1_pt_=pions.at(0)->pt();
    pion1_eta_=pions.at(0)->eta();
    pion1_phi_=pions.at(0)->phi();
    pion1_charge_=pions.at(0)->charge();
 }
 
  
  if(pions_from_meson.size()==2){
    pions_dr_ = DeltaR(pions_from_meson.at(0),pions_from_meson.at(1));
    res_mass_=(pions_from_meson.at(0)->vector()+pions_from_meson.at(1)->vector()).M();
    res_pt_=(pions_from_meson.at(0)->vector()+pions_from_meson.at(1)->vector()).pt();
    res_eta_=(pions_from_meson.at(0)->vector()+pions_from_meson.at(1)->vector()).eta();
    res_phi_=(pions_from_meson.at(0)->vector()+pions_from_meson.at(1)->vector()).phi();
  } 



  tree_->Fill();
  return 0;
  }

  void HVMGenAnalysis::PrintInfo() {}

}
