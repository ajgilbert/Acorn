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
#include "Acorn/NTupler/interface/EventInfo.h"
#include "functions.h"

#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/PFJet.h"

#include "hist_maker.h"


double hist_maker(std::vector<std::string> input_files, TString outfile_, double cs, TString Dilepton, std::string data_files, int cuts_lev){



        TFile* fOut_d = new TFile(outfile_.Data(), "RECREATE");

        int iEvt = 0; //event number


        int event=0;
        int pevent = 0;
        int nevent = 0;
        int weight;
        //double cosp;

        double luminosity = 35.9e3;

        //Bins for asymmetry
        std::vector<int> bin_p(7, 0);
        std::vector<int> bin_n(7, 0);

        std::vector<double> ybinning;

        //Define Selection
        double pT_gamma = 50;
        double eta_gamma = 2.5;
        double ET_miss = 20;
        double pT_lepton = 20;
        double eta_lepton = 2.5;
        double pT_bjet = 20;
        double eta_b = 2.5;
        double eta_j = 5.0;

        //Histograms: Dilepton

        //Make Histograms: pT
        TH1F* hPhotonPt = new TH1F("hPhotonPt", Dilepton, 100, 0, 600);
        TH1F* hLeptonPt = new TH1F("hLeptonPt", Dilepton, 100, 0, 500);
        TH1F* hBJetPt = new TH1F("hBJetPt", Dilepton, 100, 0, 500);
        TH1F* hJetPt = new TH1F("hJetPt", Dilepton, 100, 0, 500);

        //Make Histograms: nJets
        TH1F* hNBJets = new TH1F("hNBJets", Dilepton, 10, 0, 10);
        TH1F* hNJets = new TH1F("hNJets", Dilepton, 10, 0, 10);

        //Make Histograms: MET
        TH1F* hMET = new TH1F("hMET", Dilepton, 100, 0, 300);

        //Make Histograms: eta
        TH1F* hPhotonEta = new TH1F("hPhotonEta", Dilepton, 50, -6, 6);
        TH1F* hLeptonEta = new TH1F("hLeptonEta", Dilepton, 50, -6, 6);
        TH1F* hBJetEta = new TH1F("hBJetEta", Dilepton, 50, -6, 6);
        TH1F* hJetEta = new TH1F("hJetEta", Dilepton, 50, -6, 6);

        //Make Histograms: phi
        TH1F* hPhotonPhi = new TH1F("hPhotonPhi", Dilepton, 50, -6, 6);
        TH1F* hLeptonPhi = new TH1F("hLeptonPhi", Dilepton, 50, -6, 6);
        TH1F* hBJetPhi = new TH1F("hBJetPhi", Dilepton, 50, -6, 6);
        TH1F* hJetPhi = new TH1F("hJetPhi", Dilepton, 50, -6, 6);

        //Make Histograms: Delta R
        TH1F* deltaR_pl = new TH1F("deltaR_pl", Dilepton, 50, 0, 10);
        TH1F* deltaR_pj = new TH1F("deltaR_pj", Dilepton, 50, 0, 10);
        TH1F* deltaR_pb = new TH1F("deltaR_pb", Dilepton, 50, 0, 10);
        TH1F* deltaR_lj = new TH1F("deltaR_lj", Dilepton, 50, 0, 10);
        TH1F* deltaR_lb = new TH1F("deltaR_lb", Dilepton, 50, 0, 10);
        TH1F* deltaR_jb = new TH1F("deltaR_jb", Dilepton, 50, 0, 10);

	//Make histograms: Number of Leptons
	TH1F* hNLeptons = new TH1F("hNLeptons", Dilepton, 10, 0, 10);


        for (auto const& fileName : input_files) {
                TFile *f = TFile::Open(fileName.c_str(), "READ");
                TTreeReader reader("EventTree", f);
                TTreeReaderValue<std::vector<ac::GenParticle>> genParts(reader, "genParticles");
                TTreeReaderValue<std::vector<ac::Candidate>> genJets(reader, "genJets");
                TTreeReaderValue<std::vector<ac::Met>> genMet(reader, "genMet");

                TTreeReaderValue<ac::EventInfo> eventInfo(reader, "eventInfo");

                while (reader.Next()) {
                        if(iEvt % 10000 == 0) std::cout << "Processing event " << iEvt << "...\n";
                        ++iEvt; //increment

                        ac::GenParticle const* neutrino = nullptr;


                        std::vector<ac::GenParticle> bquarks;
                        std::vector<ac::GenParticle> photon_holder;
                        std::vector<ac::GenParticle> photons;
                        std::vector<ac::GenParticle> leptons;

                        ac::Candidate const* met = nullptr;


                        //Count total weights in file

                        weight = eventInfo->totalWeight();

                        if(weight >= 0){
                                pevent = pevent + weight;
                        }

                        if(weight < 0){
                                nevent = nevent + weight;
                        }

                        //Identify particles
                        for(ac::GenParticle const& part : *genParts){
                                if (part.status() == 1 && (abs(part.pdgId()) == 11 || abs(part.pdgId()) == 13)) {


                                        if (part.pt() > pT_lepton && abs(part.eta())< eta_lepton){

                                                leptons.push_back(part);

                                        }
                                }

                                //look for photons
                                if( abs(part.pdgId()) == 22 ) {

                                        photon_holder.push_back(part);

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
                        ac::GenParticle const* particle = nullptr;

                                for (ac::Candidate const& jet : *genJets) {
                                        bool isBTagged = false;
                                        bool isOverlapped = false;
                                        for (auto const& bquark : bquarks) {
                                                if (isBJet(jet, bquark)){
                                                        isBTagged = true;
                                                        break;
                                                }
                                        }

                                        if(isBTagged && jet.pt()>pT_bjet && abs(jet.eta()) <eta_b){
                                                      bjets.push_back(jet);
                                        }

                                        if(!isBTagged && jet.pt() > pT_bjet && abs(jet.eta())< eta_j && abs(jet.eta()) > 0.5 && cuts_lev >0){
                                                //if(jets.size() != 0){
                                                        //for(auto const& fjet : jets)
                                                                //if(check_overlap(jet, fjet) <= 0.4){
                                                                //      isOverlapped = true;
                                                                //}
                                                //      for(auto const& p : photon_holder){
                                                //              if(check_overlap(jet, p)<= 0.05){
                                                //                      isOverlapped = true;
                                                //              }
                                                //      }
                                                //}
                                                //if(isOverlapped == false){
                                                jets.push_back(jet);
                                                //}
                                        }
					if(!isBTagged && jet.pt()> pT_bjet && abs(jet.eta())< eta_j && cuts_lev == 0){
						jets.push_back(jet);
					}

                                }

                                //Delta R variables

                                double R_l_p;
                                double R_l_j;
                                double R_l_b;
                                double R_p_j;
                                double R_p_b;
                                double R_j_b;


                                ac::Candidate const* bjet = nullptr;
                                ac::Candidate const* jet = nullptr;
                                ac::GenParticle const* photon = nullptr;
                                ac::GenParticle const* lepton = nullptr;

                                //Sort vectors from highest pT to lowest, then select particles
                                //that pass overlap test

                                if(leptons.size() == 1 && photons.size() != 0 && bjets.size() != 0 && jets.size() != 0){
					if((cuts_lev > 0 && bjets.size() < 2 && jets.size() < 5) || (cuts_lev == 0)){

                                        	sort(leptons.begin(), leptons.end(), compare_pt_parts);
                                        	sort(photons.begin(), photons.end(), compare_pt_parts);
                                        	sort(bjets.begin(), bjets.end(), compare_pt_jets);
                                        	sort(jets.begin(), jets.end(), compare_pt_jets);

                                        	for (auto const& lep : leptons){
                                                	lepton = &lep;
                                                	for (auto const& pho : photons){
                                                        	if(check_overlap_part(pho, lepton) >= 0.4){

                                                                	R_l_p = check_overlap_part(pho,lepton);
                                                                	photon = &pho;

                                                                	for(auto const& b : bjets){
                                                                        	if(check_overlap_jet(b, lepton) >= 0.4 && check_overlap_jet(b, photon) >= 0.4 ){
	
        	                                                                        R_l_b = check_overlap_jet(b, lepton);
                	                                                                R_p_b = check_overlap_jet(b, photon);
                        	                                                        bjet = &b;
	
        	                                                                        for(auto const& j : jets){
                	                                                                        if(check_overlap_bjet(j, bjet) >= .4 && check_overlap_jet(j, photon) >= .4 && check_overlap_jet(j, lepton) >= 0.4){
													if((cuts_lev > 0 && abs(j.eta()) >= 0.5) || cuts_lev == 0){
	
        	                                                                                                R_l_j = check_overlap_bjet(j, lepton);
                	                                                                                        R_p_j = check_overlap_bjet(j, photon);
                        	                                                                                R_j_b = check_overlap_bjet(j, bjet);
	
        	                                                                                                jet = &j;
													}
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
				}



                        //implement final cuts. Require objects and MET cut
                        if(photon && lepton && bjet && jet && met->pt()>ET_miss){

                                event = weight + event;

                                hPhotonPt->Fill(photon->pt(), weight);
                                hLeptonPt->Fill(lepton->pt(), weight);
                                hBJetPt->Fill(bjet->pt(), weight);
                                hJetPt->Fill(jet->pt(), weight);

                                hPhotonEta->Fill(photon->eta(), weight);
                                hLeptonEta->Fill(lepton->eta(), weight);
                                hBJetEta->Fill(bjet->eta(), weight);
                                hJetEta->Fill(jet->eta(), weight);

                                hPhotonPhi->Fill(photon->phi(), weight);
                                hLeptonPhi->Fill(lepton->phi(), weight);
                                hBJetPhi->Fill(bjet->phi(), weight);
                                hJetPhi->Fill(jet->phi(), weight);


                                hNBJets->Fill(bjets.size(), weight);
                                hNJets->Fill(jets.size(), weight);

                                hMET->Fill(met->pt(), weight);

                                deltaR_pl -> Fill(R_l_p, weight);
                                deltaR_pj -> Fill(R_p_j, weight);
                                deltaR_pb -> Fill(R_p_b, weight);
                                deltaR_lj -> Fill(R_l_j, weight);
                                deltaR_lb -> Fill(R_l_b, weight);
                                deltaR_jb -> Fill(R_j_b, weight);

				hNLeptons -> Fill(leptons.size(), weight);

				//double l_theta = 2*atan(exp(-lepton->eta()));
				//double p_theta = 2*atan(exp(-photon->eta()));
				
				//if(l_theta < 0){ l_theta = l_theta + 2*M_PI;}
				//if(p_theta < 0){p_theta = p_theta + 2*M_PI;}

                                //double deltaPhi = cos(lepton->phi() - photon->phi());
                                //double deltaEta = cos(lepton->eta() - photon->eta());
                                //double asymm = sin(l_theta)*sin(p_theta)*deltaPhi + cos(p_theta)*cos(l_theta);
				//cout<<asymm<< "\n"<<endl;
				
				double asymm = ROOT::Math::VectorUtil::CosTheta(lepton->vector(), photon->vector());
                                int counter = floor(photon->pt()/100.);
	
                                if(asymm >= 0 && counter < 6 ){
                                        //double counter = photon->pt()/50;
                                        bin_p[counter]+= weight;
                                }

                                if(asymm < 0 && counter < 6 ){
                                        //double counter = photon->pt()/50;
                                        bin_n[counter]+=weight;
                                }
                        }

                }

                f->Close();
        }
	
		
        for(int i = 0; i < 6; i++){
                double numer = bin_p[i] - bin_n[i];
                double denom = bin_p[i] + bin_n[i];
                if(denom > 0){
                double asymmetry = numer/denom;
                cout<<"ASYMMETRY = "<<asymmetry<<"\n"<<endl;
		cout<<"Positive = "<<bin_p[i]<<"\nNegative = "<<bin_n[i]<<"\n"<<endl;
                ybinning.push_back(asymmetry);
        }
        }



        std::ofstream outFile;
        outFile.open(data_files, ios::out);
        for(int j = 0; j < ybinning.size(); j++){
             	  outFile <<ybinning[j]<<"\n"<<endl;
        }
        outFile.close();


        double tevent = pevent + nevent;
        double scale = luminosity*(cs/tevent);

        hPhotonPt->Scale(scale);
        hLeptonPt->Scale(scale);
        hBJetPt->Scale(scale);
        hJetPt->Scale(scale);

        hPhotonEta->Scale(scale);
        hLeptonEta->Scale(scale);
        hBJetEta->Scale(scale);
        hJetEta->Scale(scale);

        hPhotonPhi->Scale(scale);
        hLeptonPhi->Scale(scale);
        hBJetPhi->Scale(scale);
        hJetPhi->Scale(scale);

        hNBJets->Scale(scale);
        hNJets->Scale(scale);

        hMET->Scale(scale);
	
	hNLeptons->Scale(scale);

        deltaR_pl->Scale(scale);
        deltaR_pj->Scale(scale);
        deltaR_pb->Scale(scale);
        deltaR_lj->Scale(scale);
        deltaR_lb->Scale(scale);
        deltaR_jb->Scale(scale);

        double event_count = event*scale;

        cout<<"EVENTS IN FILE: "<< tevent<< "\n";
	cout<<"WEIGHTED EVENTS: "<< event<<"\n";
        cout<<"TOTAL: "<<event_count<< "\n";
        /////////////////////////////////////////////////////////////////////////

        cout<<"WRITING TO DISK \n"<<endl;
        fOut_d->cd();
        hPhotonPt->Write();
        hLeptonPt->Write();
        hBJetPt->Write();
        hJetPt->Write();

        hPhotonEta->Write();
        hLeptonEta->Write();
        hBJetEta->Write();
        hJetEta->Write();


        hPhotonPhi->Write();
        hLeptonPhi->Write();
        hBJetPhi->Write();
        hJetPhi->Write();

        hMET->Write();

        hNBJets->Write();
        hNJets->Write();
	hNLeptons->Write();

        deltaR_pl->Write();
        deltaR_pj->Write();
        deltaR_pb->Write();
        deltaR_lj->Write();
        deltaR_lb->Write();
        deltaR_jb->Write();


        fOut_d->Write();
        fOut_d->Close();


        return event_count;

}

double hist_maker_reco(std::vector<std::string> input_files, TString outfile_, double cs, TString Dilepton, std::string data_files, int cuts_lev){



        TFile* fOut_d = new TFile(outfile_.Data(), "RECREATE");

        int iEvt = 0; //event number


        int event=0;
        int pevent = 0;
        int nevent = 0;
        int weight;

        double luminosity = 35.9e3;

        //Bins for asymmetry
        std::vector<int> bin_p(7, 0);
        std::vector<int> bin_n(7, 0);

        std::vector<double> ybinning;

        //Define Selection
        double pT_gamma = 50;
        double eta_gamma = 2.5;
        double ET_miss = 20;
        double pT_lepton = 20;
        double eta_lepton = 2.5;
        double pT_bjet = 20;
        double eta_b = 2.5;
        double eta_j = 5.0;

        //Histograms: Dilepton

        //Make Histograms: pT
        TH1F* hPhotonPt = new TH1F("hPhotonPt", Dilepton, 150, 0, 600);
        TH1F* hLeptonPt = new TH1F("hLeptonPt", Dilepton, 150, 0, 450);
        TH1F* hBJetPt = new TH1F("hBJetPt", Dilepton, 150, 0, 450);
        TH1F* hJetPt = new TH1F("hJetPt", Dilepton, 150, 0, 450);

        //Make Histograms: nJets
        TH1F* hNBJets = new TH1F("hNBJets", Dilepton, 10, 0, 10);
        TH1F* hNJets = new TH1F("hNJets", Dilepton, 10, 0, 10);

        //Make Histograms: MET
        TH1F* hMET = new TH1F("hMET", Dilepton, 150, 0, 300);

        //Make Histograms: eta
        TH1F* hPhotonEta = new TH1F("hPhotonEta", Dilepton, 50, -6, 6);
        TH1F* hLeptonEta = new TH1F("hLeptonEta", Dilepton, 50, -6, 6);
        TH1F* hBJetEta = new TH1F("hBJetEta", Dilepton, 50, -6, 6);
        TH1F* hJetEta = new TH1F("hJetEta", Dilepton, 50, -6, 6);

        //Make Histograms: phi
        TH1F* hPhotonPhi = new TH1F("hPhotonPhi", Dilepton, 50, -6, 6);
        TH1F* hLeptonPhi = new TH1F("hLeptonPhi", Dilepton, 50, -6, 6);
        TH1F* hBJetPhi = new TH1F("hBJetPhi", Dilepton, 50, -6, 6);
        TH1F* hJetPhi = new TH1F("hJetPhi", Dilepton, 50, -6, 6);

        //Make Histograms: Delta R
        TH1F* deltaR_pl = new TH1F("deltaR_pl", Dilepton, 50, 0, 10);
        TH1F* deltaR_pj = new TH1F("deltaR_pj", Dilepton, 50, 0, 10);
        TH1F* deltaR_pb = new TH1F("deltaR_pb", Dilepton, 50, 0, 10);
        TH1F* deltaR_lj = new TH1F("deltaR_lj", Dilepton, 50, 0, 10);
        TH1F* deltaR_lb = new TH1F("deltaR_lb", Dilepton, 50, 0, 10);
        TH1F* deltaR_jb = new TH1F("deltaR_jb", Dilepton, 50, 0, 10);

	//Make Histograms: Number of Electrons
	TH1F* hNLeptons = new TH1F("hNLeptons", Dilepton, 10, 0, 10);

        for (auto const& fileName : input_files) {
                TFile *f = TFile::Open(fileName.c_str(), "READ");
                TTreeReader reader("EventTree", f);
                
		TTreeReaderValue<std::vector<ac::Photon>> photons(reader, "photons");
                TTreeReaderValue<std::vector<ac::PFJet>> pfJets(reader, "pfJets");
                TTreeReaderValue<std::vector<ac::Met>> puppiMet(reader, "puppiMet");
		TTreeReaderValue<std::vector<ac::Muon>> muons(reader, "muons");
		TTreeReaderValue<std::vector<ac::Electron>> electrons(reader, "electrons");
			
                TTreeReaderValue<ac::EventInfo> eventInfo(reader, "eventInfo");

                while (reader.Next()) {
                        if(iEvt % 10000 == 0) std::cout << "Processing event " << iEvt << "...\n";
                        ++iEvt; //increment


			//Count total weights in file

                        weight = eventInfo->totalWeight();

                        if(weight >= 0){
                                pevent = pevent + weight;
                        }

                        if(weight < 0){
                                nevent = nevent + weight;
                        }

			std::vector<ac::Electron> Electrons;
			std::vector<ac::Muon> Muons;
			std::vector<ac::Photon> Photons;
			std::vector<ac::PFJet> Jets;
			std::vector<ac::PFJet> BJets;

			//Identify particles
			 
			for( ac::Electron const& elec : *electrons){
 				double e_iso_1 = 0.0478 + 0.506/elec.pt();
				double e_iso_2 = 0.0658 + 0.963/elec.pt();

		  		if (elec.pt() > pT_lepton && abs(elec.eta())< eta_lepton  && elec.isCutBasedMediumElectron() && (fabs(elec.scEta()) < 1.4442 || fabs(elec.scEta())) ){
					if(abs(elec.eta()) <= 1.479 && elec.relativeEAIso() < e_iso_1 && fabs(elec.dxy()) < 0.05 && fabs(elec.dz()) < 0.1){
              				 	Electrons.push_back(elec);
					}
					if(abs(elec.eta())> 1.479 && elec.relativeEAIso() < e_iso_2 && fabs(elec.dxy()) < 0.10 && fabs(elec.dz()) < 0.2){
						Electrons.push_back(elec);
					}

                                 }
                   	}

			for( ac::Muon const& mu : *muons){
				double mu_iso = (mu.pfIsoSumChargedHadronPt() + max(0., mu.pfIsoSumNeutralHadronEt() + mu.pfIsoSumPhotonEt() - 0.5 + mu.pfIsoSumPUPt()))/mu.pt();
				if (mu.pt() > pT_lepton && abs(mu.eta()) < eta_lepton && mu.isMediumMuon() &&  fabs(mu.dxy())< 0.05 && fabs(mu.dz()) < 0.2){

					Muons.push_back(mu);
				}
			}

                        //look for photons
			for( ac::Photon const& pho : *photons){
				if(pho.pt() > pT_gamma && abs(pho.eta()) < eta_gamma && pho.isMediumIdPhoton() && !pho.hasPixelSeed() && (fabs(pho.scEta()) < 1.4442 || fabs(pho.scEta()) > 1.566)){

					Photons.push_back(pho);
				}
			}
			
			for( ac::PFJet const& jet_i : *pfJets){
				if(jet_i.deepCSVDiscriminatorBvsAll() > 0.7527 && jet_i.pt() > pT_bjet && abs(jet_i.eta()) < eta_b && jet_i.passesJetID()){
					BJets.push_back(jet_i);
				}
				if(jet_i.deepCSVDiscriminatorBvsAll() <= 0.7527 && jet_i.pt() > pT_bjet && abs(jet_i.eta()) < eta_j && jet_i.passesJetID()){
					Jets.push_back(jet_i);
				}
			}


                        auto met = puppiMet->begin();


                                //Delta R variables

                                double R_e_p;
                                double R_e_j;
                                double R_e_b;
                                double R_p_j;
                                double R_p_b;
                                double R_j_b;
				
				double R_m_p_2;
                                double R_m_j_2;
                                double R_m_b_2;
                                double R_p_j_2;
                                double R_p_b_2;
                                double R_j_b_2;

                                ac::PFJet const* bjet = nullptr;
                                ac::PFJet const* jet = nullptr;
                                ac::Photon const* photon = nullptr;
                                ac::Electron const* electron = nullptr;

                                //Sort vectors from highest pT to lowest, then select particles
                                //that pass overlap test

                                if((Electrons.size() + Muons.size() == 1) && Photons.size() != 0 && BJets.size() != 0 && Jets.size() != 0){
                                        if((cuts_lev >= 1 && BJets.size() < 2 && Jets.size() < 5) || (cuts_lev == 0)){

						sort(Electrons.begin(), Electrons.end(), compare_pt_e);
						sort(Photons.begin(), Photons.end(), compare_pt_p);
                                                sort(BJets.begin(), BJets.end(), compare_pt_j);
                                                sort(Jets.begin(), Jets.end(), compare_pt_j);

                                                for (auto const& ele : Electrons){
                                                        electron = &ele;
                                                        for (auto const& pho : Photons){
                                                                if(check_overlap_pe(pho, electron) >= 0.4){

                                                                        R_e_p = check_overlap_pe(pho,electron);
                                                                        photon = &pho;

                                                                        for(auto const& b : BJets){
                                                                                if(check_overlap_be(b, electron) >= 0.4 && check_overlap_bp(b, photon) >= 0.4 ){

                                                                                        R_e_b = check_overlap_be(b, electron);
                                                                                        R_p_b = check_overlap_bp(b, photon);
                                                                                        bjet = &b;

                                                                                        for(auto const& j : Jets){
                                                                                                if(check_overlap_jb(j, bjet) >= .4 && check_overlap_jp(j, photon) >= .4 && check_overlap_je(j, electron) >= 0.4){// && abs(j.eta()) >= 0.5){
													if((cuts_lev > 0 && abs(j.eta()) >= 0.5) || cuts_lev == 0){

                                                                                                                R_e_j = check_overlap_je(j, electron);
                                                                                                                R_p_j = check_overlap_jp(j, photon);
                                                                                                                R_j_b = check_overlap_jb(j, bjet);

                                                                                                                jet = &j;
													}
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
                                }


			ac::PFJet const* bjet_2 = nullptr;
                        ac::PFJet const* jet_2 = nullptr;
                        ac::Photon const* photon_2 = nullptr;
                        ac::Muon const* muon_2 = nullptr;


			if((Electrons.size() + Muons.size() == 1) && Photons.size() != 0 && BJets.size() != 0 && Jets.size() != 0){
                                        if((cuts_lev > 0 && BJets.size() < 2 && Jets.size() < 5) || (cuts_lev == 0)){

                                                sort(Muons.begin(), Muons.end(), compare_pt_m);
                                                sort(Photons.begin(), Photons.end(), compare_pt_p);
                                                sort(BJets.begin(), BJets.end(), compare_pt_j);
                                                sort(Jets.begin(), Jets.end(), compare_pt_j);

                                                for (auto const& mu : Muons){
                                                        muon_2 = &mu;
                                                        for (auto const& pho : Photons){
                                                                if(check_overlap_pm(pho, muon_2) >= 0.4){

                                                                        R_m_p_2 = check_overlap_pm(pho, muon_2);
                                                                        photon_2 = &pho;

                                                                        for(auto const& b : BJets){
                                                                                if(check_overlap_bm(b, muon_2) >= 0.4 && check_overlap_bp(b, photon_2) >= 0.4 ){

                                                                                        R_m_b_2 = check_overlap_bm(b, muon_2);
                                                                                        R_p_b_2 = check_overlap_bp(b, photon_2);
                                                                                        bjet_2 = &b;

                                                                                        for(auto const& j : Jets){
                                                                                                if(check_overlap_jb(j, bjet_2) >= .4 && check_overlap_jp(j, photon_2) >= .4 && check_overlap_jm(j, muon_2) >= 0.4){// && abs(j.eta()) >= 0.5){
													if((cuts_lev > 0 && abs(j.eta()) >= 0.5) || cuts_lev == 0){

                                                                                                                R_m_j_2 = check_overlap_jm(j, muon_2);
                                                                                                                R_p_j_2 = check_overlap_jp(j, photon_2);
                                                                                                                R_j_b_2 = check_overlap_jb(j, bjet_2);

                                                                                                                jet_2 = &j;
                                                                                                
													}
												}
                                                                                                if(jet_2){break;}
                                                                                        }
                                                                                }
                                                                                if(jet_2){break;}
                                                                        }
                                                                }
                                                                if(jet_2){break;}
                                                        }
                                                        if(jet_2){break;}
                                                }
                                        }
                                }



			int muEvent = 0;
			int elEvent = 0;
			
                        //implement final cuts. Require objects and MET cut
                        if(photon && electron && bjet && jet && met->pt()>ET_miss && (!jet_2)){
				if(cuts_lev ==0 || (cuts_lev >0 && abs(jet->eta()) >= 0.5)){
					elEvent = 1;
				}
			}
			if(photon_2 && muon_2 && bjet_2 && jet_2 && met->pt()>ET_miss && (!jet)){
				if(cuts_lev ==0 || (cuts_lev > 0 && abs(jet_2->eta()) >= 0.5)){
					muEvent = 1;
				}
			}
			if(photon && electron && bjet && jet && photon_2 && muon_2 && bjet_2 && jet_2 && met->pt()>ET_miss){
				if( electron->pt() > muon_2->pt()){
					
					int elEvent = 1;
			}
				else{
					int muEvent = 1;
				}
			}
			
			if(elEvent){
                                event = weight + event;

                                hPhotonPt->Fill(photon->pt(), weight);
                                hLeptonPt->Fill(electron->pt(), weight);
                                hBJetPt->Fill(bjet->pt(), weight);
                                hJetPt->Fill(jet->pt(), weight);

                                hPhotonEta->Fill(photon->eta(), weight);
                                hLeptonEta->Fill(electron->eta(), weight);
                                hBJetEta->Fill(bjet->eta(), weight);
                                hJetEta->Fill(jet->eta(), weight);

                                hPhotonPhi->Fill(photon->phi(), weight);
                                hLeptonPhi->Fill(electron->phi(), weight);
                                hBJetPhi->Fill(bjet->phi(), weight);
                                hJetPhi->Fill(jet->phi(), weight);


                                hNBJets->Fill(BJets.size(), weight);
                                hNJets->Fill(Jets.size(), weight);

                                hMET->Fill(met->pt(), weight);

				hNLeptons-> Fill(Electrons.size()+Muons.size(), weight);

                                deltaR_pl -> Fill(R_e_p, weight);
                                deltaR_pj -> Fill(R_p_j, weight);
                                deltaR_pb -> Fill(R_p_b, weight);
                                deltaR_lj -> Fill(R_e_j, weight);
                                deltaR_lb -> Fill(R_e_b, weight);
				deltaR_jb -> Fill(R_j_b, weight);

				//double l_theta = 2*atan(exp(-electron->eta()));
                                //double p_theta = 2*atan(exp(-photon->eta()));


                                //double deltaPhi = cos(electron->phi() - photon->phi());
                                //double asymm = sin(l_theta)*sin(p_theta)*deltaPhi + cos(p_theta)*cos(l_theta);
				
				double asymm = ROOT::Math::VectorUtil::CosTheta(electron->vector(), photon->vector());
                                int counter = floor(photon->pt()/100.);

                                if(asymm >= 0 && counter < 6){
                                        //double counter = photon->pt()/50;
                                        bin_p[counter]++;
                                }

                                if(asymm < 0 && counter < 6){
                                        //double counter = photon->pt()/50;
                                        bin_n[counter]++;
                                }
                        }

			if(muEvent){
                                event = weight + event;

                                hPhotonPt->Fill(photon_2->pt(), weight);
                                hLeptonPt->Fill(muon_2->pt(), weight);
                                hBJetPt->Fill(bjet_2->pt(), weight);
                                hJetPt->Fill(jet_2->pt(), weight);

                                hPhotonEta->Fill(photon_2->eta(), weight);
                                hLeptonEta->Fill(muon_2->eta(), weight);
                                hBJetEta->Fill(bjet_2->eta(), weight);
                                hJetEta->Fill(jet_2->eta(), weight);

                                hPhotonPhi->Fill(photon_2->phi(), weight);
                                hLeptonPhi->Fill(muon_2->phi(), weight);
                                hBJetPhi->Fill(bjet_2->phi(), weight);
                                hJetPhi->Fill(jet_2->phi(), weight);


                                hNBJets->Fill(BJets.size(), weight);
                                hNJets->Fill(Jets.size(), weight);
				hNLeptons->Fill(Electrons.size()+Muons.size(), weight);

                                hMET->Fill(met->pt(), weight);
				
                                deltaR_pl -> Fill(R_m_p_2, weight);
                                deltaR_pj -> Fill(R_p_j_2, weight);
                                deltaR_pb -> Fill(R_p_b_2, weight);
                                deltaR_lj -> Fill(R_m_j_2, weight);
                                deltaR_lb -> Fill(R_m_b_2, weight);
                                deltaR_jb -> Fill(R_j_b_2, weight);

				//double l_theta = 2*atan(exp(-muon_2->eta()));
                                //double p_theta = 2*atan(exp(-photon_2->eta()));

                                //if(l_theta < 0){ l_theta = l_theta + 2*M_PI;}
                                //if(p_theta < 0){p_theta = p_theta + 2*M_PI;}

                                //double deltaPhi = cos(muon_2->phi() - photon_2->phi());
                                //double deltaEta = cos(lepton->eta() - photon->eta());
                                //double asymm = sin(l_theta)*sin(p_theta)*deltaPhi + cos(p_theta)*cos(l_theta);

                                //if(asymm > 1){
                                //        cout<<"ASYMMETRY GREATER THAN ONE ------------------------------------------- \n"<<endl;
                                //        return 0;
                                //}

				double asymm = ROOT::Math::VectorUtil::CosTheta(muon_2->vector(), photon_2->vector());

                                int counter = floor(photon_2->pt()/100.);

                                if(asymm >= 0 && counter < 6){
                                        //double counter = photon->pt()/50;
                                        bin_p[counter]++;
                                }

                                if(asymm < 0 && counter < 6){
                                        //double counter = photon->pt()/50;
                                        bin_n[counter]++;
                                }
                        }




                }

                f->Close();
        }


        for(int i = 0; i < 6; i++){
                double numer = bin_p[i] - bin_n[i];
                double denom = bin_p[i] + bin_n[i];
                if(denom > 0){
                double asymmetry = numer/denom;
		cout<<"ASYMMETRY = "<< asymmetry <<"\n"<< endl;
                cout<<"Positive = "<<bin_p[i]<<"\nNegative = "<<bin_n[i]<<"\n"<<endl;
                ybinning.push_back(asymmetry);
        }
        }



        std::ofstream outFile;
        outFile.open(data_files, ios::out);
        for(int j = 0; j < ybinning.size(); j++){
                outFile <<ybinning[j]<<"\n"<<endl;
        }
        outFile.close();


        double tevent = pevent + nevent;
        double scale = luminosity*(cs/tevent);

        hPhotonPt->Scale(scale);
        hLeptonPt->Scale(scale);
        hBJetPt->Scale(scale);
        hJetPt->Scale(scale);

        hPhotonEta->Scale(scale);
        hLeptonEta->Scale(scale);
        hBJetEta->Scale(scale);
        hJetEta->Scale(scale);

        hPhotonPhi->Scale(scale);
        hLeptonPhi->Scale(scale);
        hBJetPhi->Scale(scale);
        hJetPhi->Scale(scale);

        hNBJets->Scale(scale);
        hNJets->Scale(scale);

	hNLeptons->Scale(scale);

        hMET->Scale(scale);

        deltaR_pl->Scale(scale);
        deltaR_pj->Scale(scale);
        deltaR_pb->Scale(scale);
        deltaR_lj->Scale(scale);
        deltaR_lb->Scale(scale);
        deltaR_jb->Scale(scale);

        double event_count = event*scale;

        cout<<"EVENTS IN FILE: "<< tevent<< "\n";
        cout<<"TOTAL: "<<event_count<< "\n";
        /////////////////////////////////////////////////////////////////////////

        cout<<"WRITING TO DISK \n"<<endl;
        fOut_d->cd();
        hPhotonPt->Write();
        hLeptonPt->Write();
        hBJetPt->Write();
        hJetPt->Write();

        hPhotonEta->Write();
        hLeptonEta->Write();
        hBJetEta->Write();
        hJetEta->Write();


        hPhotonPhi->Write();
        hLeptonPhi->Write();
        hBJetPhi->Write();
        hJetPhi->Write();
	hJetPhi->Write();

        hMET->Write();

        hNBJets->Write();
        hNJets->Write();
	hNLeptons->Write();

        deltaR_pl->Write();
        deltaR_pj->Write();
        deltaR_pb->Write();
        deltaR_lj->Write();
        deltaR_lb->Write();
        deltaR_jb->Write();


        fOut_d->Write();
        fOut_d->Close();


        return event_count;

}

