#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <vector>
#include <string>
#include <iostream>
#include <math.h>

#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/Met.h"

//TString outfile_ = "basic_analysis_3.hist";

bool isBJet(ac::Candidate const& jet, ac::GenParticle const& bquark){
	double delta_R_max = .3;
	TLorentzVector jet_lv, b_lv;
	jet_lv.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.M());
	b_lv.SetPtEtaPhiM(bquark.pt(), bquark.eta(), bquark.phi(), bquark.M());
	bool retval = abs(jet_lv.DeltaR(b_lv))< delta_R_max;
	//std::cout<<"isBjet = " << retval << std::endl;
	return retval;
}

//functions to check overlap.  If overlap -> return true

bool check_overlap_part(ac::GenParticle const& part_1, ac::GenParticle const* part_2){
        double delta_R_max = .4;
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        bool overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))) < delta_R_max;
        return overlap;
}

bool check_overlap_jet(ac::Candidate const& jet_1, ac::GenParticle const* part_2){
        double delta_R_max = .4;
        TLorentzVector jet_lv_1, part_lv_2;
        jet_lv_1.SetPtEtaPhiM(jet_1.pt(), jet_1.eta(), jet_1.phi(), jet_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        bool overlap = abs(sqrt((jet_lv_1.Phi() - part_lv_2.Phi())*(jet_lv_1.Phi() - part_lv_2.Phi()) + (jet_lv_1.Eta() - part_lv_2.Eta())*(jet_lv_1.Eta() - part_lv_2.Eta()))) < delta_R_max;
        return overlap;
}

bool check_overlap_bjet(ac::Candidate const& jet_1, ac::Candidate const* jet_2){
        double delta_R_max = .4;
        TLorentzVector jet_lv_1, jet_lv_2;
        jet_lv_1.SetPtEtaPhiM(jet_1.pt(), jet_1.eta(), jet_1.phi(), jet_1.M());
        jet_lv_2.SetPtEtaPhiM(jet_2->pt(), jet_2->eta(), jet_2->phi(), jet_2->M());
        bool overlap = abs(sqrt((jet_lv_1.Phi() - jet_lv_2.Phi())*(jet_lv_1.Phi() - jet_lv_2.Phi()) + (jet_lv_1.Eta() - jet_lv_2.Eta())*(jet_lv_1.Eta() - jet_lv_2.Eta()))) < delta_R_max;
        return overlap;
}


//Functions to select highest pT particle; Not used


ac::GenParticle const* select_highest_particle(std::vector<ac::GenParticle> particles){
        ac::GenParticle const* particle = nullptr;

        for (auto const& p : particles)  {
                                if(!particle || particle->pt() < p.pt()){

                                                particle = &p;

                                }
                        }
        return particle;
}


ac::GenParticle const* select_particle(std::vector<ac::GenParticle> parts, ac::GenParticle const* ini_part){
        ac::GenParticle const* particle = nullptr;
	
	if(!ini_part){
        	for (auto const& p : parts)  {
								
			if(!particle||(particle->pt() < p.pt())){
			
				particle = &p;
			}
		}
	}
	else{
		for (auto const& p : parts)  {

                        if(!particle||(particle->pt() < p.pt())){
			
				if(p.pt() <= ini_part->pt()){ 	
                                        particle = &p;
				}
                        }
		}
		
	}
        return particle;
}

ac::Candidate const* select_jet(std::vector<ac::Candidate> jets, ac::Candidate const* ini_jet){
        ac::Candidate const* jet = nullptr;

       	for (auto const& j : jets)  {

		if(!jet || jet->pt()<j.pt()){
			
			if(!ini_jet){
				jet = &j;
			}
			else{
				if(j.pt()<=ini_jet->pt()){
					jet = &j;
				}		
			}
		}
	}
	
        return jet;
}


//Functions used for sorting by descending pT
bool compare_pt_parts(ac::GenParticle const& part_1, ac::GenParticle const& part_2){
	TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

	return (part_lv_1.Pt()> part_lv_2.Pt());
}	

bool compare_pt_jets(ac::Candidate const& part_1, ac::Candidate const& part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

        return (part_lv_1.Pt()> part_lv_2.Pt());
}







////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////MAIN ANALYSIS FUNCTION//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int initial_distributions(){
	std::vector<std::string> input_files = {"/home/mmackenz/storage/TGJets/EventTree_1.root"};

	int iEvt = 0; //event number
	
	//Define Selection
	double pT_gamma = 50;
	double eta_gamma = 2.5;
	double ET_miss = 20;
	double pT_lepton = 20;
	double eta_lepton = 2.5;
	double pT_bjet = 20;
	double eta_b = 2.5;
	double eta_j = 5.0;	


	//Make Histograms: pT
	TH1F* hPhotonPt = new TH1F("hPhotonPt", "Leading Photon pT", 150, 0, 300);
	TH1F* hLeptonPt = new TH1F("hLeptonPt", "Leading Lepton pT", 150, 0, 300);
	TH1F* hBJetPt = new TH1F("hBJetPt", "Leading b-tagged Jet pT", 150, 0, 300);
	TH1F* hJetPt = new TH1F("hJetPt", "Leading Jet pT", 150, 0, 300);

	//Make Histograms: nJets
	TH1F* hNBJets = new TH1F("hNBJets", "N(b-tagged Jets)", 10, 0, 10);
	TH1F* hNJets = new TH1F("hNJets", "N(untagged Jets)", 10, 0, 10);

	//Make Histograms: MET
	TH1F* hMET = new TH1F("hMET", "MET", 150, 0, 300);

	//Make Histograms: eta
	TH1F* hPhotonEta = new TH1F("hPhotonEta", "Leading Photon Eta", 50, -6, 6);
	TH1F* hLeptonEta = new TH1F("hLeptonEta", "Leading Lepton Eta", 50, -6, 6);
	TH1F* hBJetEta = new TH1F("hBJetEta", "Leading b-tagged Jet Eta", 50, -6, 6);
	TH1F* hJetEta = new TH1F("hJetEta", "Leading Jet Eta", 50, -6, 6);

	//Make Histograms: phi
	TH1F* hPhotonPhi = new TH1F("hPhotonPhi", "Leading Photon Phi", 50, -6, 6);
        TH1F* hLeptonPhi = new TH1F("hLeptonPhi", "Leading Lepton Phi", 50, -6, 6);
        TH1F* hBJetPhi = new TH1F("hBJetPhi", "Leading b-tagged Jet Phi", 50, -6, 6);
        TH1F* hJetPhi = new TH1F("hJetPhi", "Leading Jet Phi", 50, -6, 6);

	

	//Histograms with Selection
	//Pt
	//TH1F* hPhotonPt_s = new TH1F("hPhotonPt_s", "Leading Photon pT (Selection)", 150, 0, 300);
        //TH1F* hLeptonPt_s = new TH1F("hLeptonPt_s", "Leading Lepton pT (Selection)", 150, 0, 300);
        //TH1F* hBJetPt_s = new TH1F("hBJetPt_s", "Leading b-tagged Jet pT (Selection)", 150, 0, 300);
        //TH1F* hJetPt_s = new TH1F("hJetPt_s", "Leading Jet pT (Selection)", 150, 0, 300);

	//nJets
        //TH1F* hNBJets_s = new TH1F("hNBJets_s", "N(b-tagged Jets) (Selection)", 10, 0, 10);
        //TH1F* hNJets_s = new TH1F("hNJets_s", "N(untagged Jets) (Selection)", 10, 0, 10);

	//MET
        //TH1F* hMET_s = new TH1F("hMET_s", "MET (Selection)", 150, 0, 300);

	//eta
	//TH1F* hPhotonEta_s = new TH1F("hPhotonEta_s", "Leading Photon Eta (Selection)", 50, -6, 6);
        //TH1F* hLeptonEta_s = new TH1F("hLeptonEta_s", "Leading Lepton Eta (Selection)", 50, -6, 6);
        //TH1F* hBJetEta_s = new TH1F("hBJetEta_s", "Leading b-tagged Jet Eta (Selection)", 50, -6, 6);
        //TH1F* hJetEta_s = new TH1F("hJetEta_s", "Leading Jet Eta (Selection)", 50, -6, 6);


	//phi
	//TH1F* hPhotonPhi_s = new TH1F("hPhotonPhi_s", "Leading Photon Phi (Selection)", 50, -6, 6);
        //TH1F* hLeptonPhi_s = new TH1F("hLeptonPhi_s", "Leading Lepton Phi (Selection)", 50, -6, 6);
        //TH1F* hBJetPhi_s = new TH1F("hBJetPhi_s", "Leading b-tagged Jet Phi (Selection)", 50, -6, 6);
        //TH1F* hJetPhi_s = new TH1F("hJetPhi_s", "Leading Jet Phi (Selection)", 50, -6, 6);




	for (auto const& fileName : input_files) {
		TFile *f = TFile::Open(fileName.c_str(), "READ");
		TTreeReader reader("EventTree", f);
		TTreeReaderValue<std::vector<ac::GenParticle>> genParts(reader, "genParticles");
		TTreeReaderValue<std::vector<ac::Candidate>> genJets(reader, "genJets");
		TTreeReaderValue<std::vector<ac::Met>> genMet(reader, "genMet");
		
		
		while (reader.Next()) {
			if(iEvt % 10000 == 0) std::cout << "Processing event " << iEvt << "...\n";
			++iEvt; //increment
			
			ac::GenParticle const* neutrino = nullptr;

			
			std::vector<ac::GenParticle> bquarks;
			
			std::vector<ac::GenParticle> photons;
                        std::vector<ac::GenParticle> leptons;
	
			ac::Candidate const* met = nullptr;
			
			//for (auto const& part : *genParts){
			for(ac::GenParticle const& part : *genParts){
				//part.Print();
				//look for muons and electrons
				if (part.status() == 1 && (abs(part.pdgId()) == 11 || abs(part.pdgId()) == 13)) {
				

					if (part.pt() > pT_lepton && abs(part.eta())< eta_lepton){ 	
						
						leptons.push_back(part);
						
					}
				}

				//look for photons
                                if( abs(part.pdgId()) == 22 ) {
					
                                        if(part.pt() > pT_gamma && abs(part.eta())< eta_gamma){
                                                
						photons.push_back(part);
						
                                     	if(iEvt == 4419){part.Print();} 
					}
                                }

				if (part.status() ==1 && (abs(part.pdgId()) == 12 || abs(part.pdgId()) == 14 || abs(part.pdgId()) == 16)) {
						
					if (!neutrino || (part.pt() > neutrino -> pt())){

                                                neutrino = &part;}

					

				
				}
				//look for b quarks
				if(abs(part.pdgId()) == 5) {
					bquarks.push_back(part);
				}
			
			}


			
			for(auto const& m : *genMet){
				met = &m;
			}
				

			//end gen particle loop
			//begin b-jet identification
			std::vector<ac::Candidate> bjets;
			std::vector<ac::Candidate> jets;
				for (ac::Candidate const& jet : *genJets) {
					bool isBTagged = false;
					for (auto const& bquark : bquarks) {
						if (isBJet(jet, bquark)){
							isBTagged = true;
							break;
						}
					}
					if(jet.pt()> pT_bjet && abs(jet.eta())< eta_j){
						jets.push_back(jet);}
					if(isBTagged && jet.pt()>pT_bjet && abs(jet.eta()) <eta_b){
						bjets.push_back(jet);
					}
					
				}
				


				//ac::GenParticle const* temp_lepton = nullptr;
				//ac::GenParticle const* temp_photon = nullptr;
				//ac::Candidate const* temp_bjet = nullptr;
				//ac::Candidate const* temp_jet = nullptr;
				
				
				
				ac::Candidate const* bjet = nullptr;
	                        ac::Candidate const* jet = nullptr;
        	                ac::GenParticle const* photon = nullptr;
                       		ac::GenParticle const* lepton = nullptr;
		
				//Sort vectors from highest pT to lowest, then select particles
				//that pass overlap test

				if(leptons.size() != 0 && photons.size() != 0 && bjets.size() != 0 && jets.size() != 0){
					sort(leptons.begin(), leptons.end(), compare_pt_parts);
                                        sort(photons.begin(), photons.end(), compare_pt_parts);
                                        sort(bjets.begin(), bjets.end(), compare_pt_jets);
                                        sort(jets.begin(), jets.end(), compare_pt_jets);

					for (auto const& lep : leptons){
						lepton = &lep;
						for (auto const& pho : photons){
							if(!check_overlap_part(pho, lepton)){
								photon = &pho;
								for(auto const& b : bjets){
									if(!check_overlap_jet(b, lepton) && !check_overlap_jet(b, photon)){
										bjet = &b;
										for(auto const& j : jets){
											if(!check_overlap_bjet(j, bjet) && !check_overlap_jet(j, photon) && !check_overlap_jet(j, lepton)){
													jet = &j;
											}
											if(jet){break;}
										}
									}
									if(jet){break;}
								}
							}
							if(jet){break;}
						}
						if(jet){break;}
					}
				}
	
			


		
						
			

			//implement final cuts. Require objects and MET cut
			if(photon && lepton && bjet && jet && met->pt()>ET_miss){			
			

				hPhotonPt->Fill(photon->pt());
				hLeptonPt->Fill(lepton->pt());
				hBJetPt->Fill(bjet->pt());
				hJetPt->Fill(jet->pt());

				hPhotonEta->Fill(photon->eta());
				hLeptonEta->Fill(lepton->eta());
				hBJetEta->Fill(bjet->eta());
				hJetEta->Fill(jet->eta());

				hPhotonPhi->Fill(photon->phi());
                                hLeptonPhi->Fill(lepton->phi());
                                hBJetPhi->Fill(bjet->phi());
                                hJetPhi->Fill(jet->phi());


				hNBJets->Fill(bjets.size());
				hNJets->Fill(jets.size());

				hMET->Fill(met->pt());
				       
		
			
			
			}
		}
			

		

		
		f->Close();
	}


	//Draw histograms
	
	//TCanvas* c8 = new TCanvas("c8", "c8", 1000,1000);
	//hLeptonPt_s->Draw();
	//hLeptonPt_s->GetXaxis()->SetTitle("pT (GeV)");

	TCanvas* c9 = new TCanvas("c9", "c9", 1000,1000);
	hPhotonPt->Draw();
	hPhotonPt->GetXaxis()->SetTitle("pT (GeV)");

	//TCanvas* c10 = new TCanvas("c10", "c10", 1000,1000);
	//hPhotonPt_s->Draw();
	//hPhotonPt_s->GetXaxis()->SetTitle("pT (GeV)");

	TCanvas* c11 = new TCanvas("c11", "c11", 1000,1000);
	hLeptonPt->Draw();
	hLeptonPt->GetXaxis()->SetTitle("pT (GeV)");

	//TCanvas* c12 = new TCanvas("c12", "c12", 1000,1000);
	//hLeptonPt_s->Draw();
	//hLeptonPt_s->GetXaxis()->SetTitle("pT (GeV)");
        
	TCanvas* c13 = new TCanvas("c13", "c13", 1000,1000);
	hBJetPt->Draw();
	hBJetPt->GetXaxis()->SetTitle("pT (GeV)");
        
	//TCanvas* c14 = new TCanvas("c14", "c14", 1000,1000);
	//hBJetPt_s->Draw();
	//hBJetPt_s->GetXaxis()->SetTitle("pT (GeV)");
	
	TCanvas* c15 = new TCanvas("c15", "c15", 1000,1000);
	hJetPt->Draw();
	hJetPt->GetXaxis()->SetTitle("pT (GeV)");
        
	//TCanvas* c16 = new TCanvas("c16", "c16", 1000,1000);
	//hJetPt_s->Draw();
	//hJetPt_s->GetXaxis()->SetTitle("pT (GeV)");
        
	TCanvas* c17 = new TCanvas("c17", "c17", 1000,1000);
	hPhotonEta->Draw();
	hPhotonEta->GetXaxis()->SetTitle("eta");
	
	//TCanvas* c18 = new TCanvas("c18", "c18", 1000,1000);
	//hPhotonEta_s->Draw();
	//hPhotonEta_s->GetXaxis()->SetTitle("eta");
        
	TCanvas* c19 = new TCanvas("c19", "c19", 1000,1000);
	hLeptonEta->Draw();
	hLeptonEta->GetXaxis()->SetTitle("eta");
        
	//TCanvas* c20 = new TCanvas("c20", "c20", 1000,1000);
	//hLeptonEta_s->Draw();
	//hLeptonEta_s->GetXaxis()->SetTitle("eta");
	
	TCanvas* c21 = new TCanvas("c21", "c21", 1000,1000);
	hBJetEta->Draw();
	hBJetEta->GetXaxis()->SetTitle("eta");
        
	//TCanvas* c22 = new TCanvas("c22", "c22", 1000,1000);
	//hBJetEta_s->Draw();
	//hBJetEta_s->GetXaxis()->SetTitle("eta");
        
	TCanvas* c23 = new TCanvas("c23", "c23", 1000,1000);
	hJetEta->Draw();
	hJetEta->GetXaxis()->SetTitle("eta");
	
	//TCanvas* c24 = new TCanvas("c24", "c24", 1000,1000);
        //hJetEta_s->Draw();
	//hJetEta_s->GetXaxis()->SetTitle("eta");

	//TCanvas* c25 = new TCanvas("c25", "c25", 1000,1000);
        //hPhotonPhi->Draw();
	//hPhotonPhi->GetXaxis()->SetTitle("phi");

        //TCanvas* c26 = new TCanvas("c26", "c26", 1000,1000);
        //hPhotonPhi_s->Draw();
	//hPhotonPhi_s->GetXaxis()->SetTitle("phi");

        //TCanvas* c27 = new TCanvas("c27", "c27", 1000,1000);
        //hLeptonPhi->Draw();
	//hLeptonPhi->GetXaxis()->SetTitle("phi");

        //TCanvas* c28 = new TCanvas("c28", "c28", 1000,1000);
        //hLeptonPhi_s->Draw();
	//hLeptonPhi_s->GetXaxis()->SetTitle("phi");

        //TCanvas* c29 = new TCanvas("c29", "c29", 1000,1000);
        //hBJetPhi->Draw();
	//hBJetPhi->GetXaxis()->SetTitle("phi");

        //TCanvas* c30 = new TCanvas("c30", "c30", 1000,1000);
        //hBJetPhi_s->Draw();
	//hBJetPhi_s->GetXaxis()->SetTitle("phi");

        //TCanvas* c31 = new TCanvas("c31", "c31", 1000,1000);
        //hJetPhi->Draw();
	//hJetPhi->GetXaxis()->SetTitle("phi");

        //TCanvas* c32 = new TCanvas("c32", "c32", 1000,1000);
        //hJetPhi_s->Draw();
	//hJetPhi_s->GetXaxis()->SetTitle("phi");

	TCanvas* c33 = new TCanvas("c33", "c33", 1000, 1000);
	hMET->Draw();
	hMET->GetXaxis()->SetTitle("MET (GeV)");

	//TCanvas* c34 = new TCanvas("c34", "c34", 1000, 1000);
	//hMET_s->Draw();
	//hMET_s->GetXaxis()->SetTitle("MET (GeV)");

	//fOut->cd();
			
	//hPhotonPt->Draw();
	//hLeptonPt->Draw();

	//fOut->Write();
	//fOut->Close();
				            
	return 0;
			

}


