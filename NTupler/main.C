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
//#include "functions.h"
#include "hist_maker.h"

int main() {
	std::vector<std::string> input_files =   {"/home/mmackenz/storage/TGJets/EventTree_1.root",//};//,
							"/home/mmackenz/storage/TGJets/EventTree_2.root",
							"/home/mmackenz/storage/TGJets/EventTree_3.root",
							"/home/mmackenz/storage/TGJets/EventTree_4.root",
							"/home/mmackenz/storage/TGJets/EventTree_5.root"};

	std::vector<std::string> input_files_d = {"/home/mmackenz/storage/TTGJets_Dilepon/EventTree_1.root",//};//,
							"/home/mmackenz/storage/TTGJets_Dilepon/EventTree_2.root",
							"/home/mmackenz/storage/TTGJets_Dilepon/EventTree_3.root",
							"/home/mmackenz/storage/TTGJets_Dilepon/EventTree_4.root",
							"/home/mmackenz/storage/TTGJets_Dilepon/EventTree_5.root"};
	
	std::vector<std::string> input_files_s = {"/home/mmackenz/storage/TTGJets_SingleLep/EventTree_1.root",//};//,
							"/home/mmackenz/storage/TTGJets_SingleLep/EventTree_2.root",
							"/home/mmackenz/storage/TTGJets_SingleLep/EventTree_3.root",
							"/home/mmackenz/storage/TTGJets_SingleLep/EventTree_4.root",
							"/home/mmackenz/storage/TTGJets_SingleLep/EventTree_5.root"};

	std::vector<TString> outfiles = {"signal_distributions_reco.hist", "dilepton_distributions_reco.hist", "singlelep_distributions_reco.hist", "signal_distributions_gen.hist", "dilepton_distributions_gen.hist", "singlelep_distributions_gen.hist"};
	std::vector<std::string> data_files = {"asymm_signal.txt", "asymm_dilepton.txt", "asymm_singlelep.txt","asymm_signal_reco.txt", "asymm_dilepton_reco.txt", "asymm_singlelep_reco.txt"};

	TString signal = "Signal";
	TString dilepton = "Dilepton";
	TString singlelep = "Single Lepton";

	double dummy4 = hist_maker(input_files, outfiles[3], 1.018, signal, data_files[0], 0);
        double dummy5 = hist_maker(input_files_d, outfiles[4], 2.3023, dilepton, data_files[1], 0);
        double dummy6 = hist_maker(input_files_s, outfiles[5], 10.50732, singlelep, data_files[2], 0);

	double dummy1 = hist_maker_reco(input_files, outfiles[0], 1.018, signal, data_files[3], 0);
	double dummy2 = hist_maker_reco(input_files_d, outfiles[1], 2.3023, dilepton, data_files[4], 0);
	double dummy3 = hist_maker_reco(input_files_s, outfiles[2], 10.50732, singlelep, data_files[5], 0);
		
	
	double sig_efficiency = dummy1/dummy4;
	double dil_efficiency = dummy2/dummy5;
	double sin_efficiency = dummy3/dummy6;
	double efficiency = (dummy1 + dummy2 + dummy3)/(dummy4 + dummy5 + dummy6);

	cout<<"signal efficiency = "<<sig_efficiency<<"\n";
	cout<<"dilepton efficiency = "<<dil_efficiency<<"\n";
	cout<<"single lepton efficiency = "<<sin_efficiency<<"\n";
	cout<<"total efficiency = "<<efficiency<<"\n";	

        //STACKS
        THStack* lepton_pt = new THStack("lepton_pt","");
        THStack* photon_pt = new THStack("photon_pt","");
        THStack* jet_pt = new THStack("jet_pt","");
        THStack* bjet_pt = new THStack("bjet_pt","");
        THStack* lepton_eta = new THStack("lepton_eta","");
        THStack* photon_eta = new THStack("photon_eta","");
        THStack* jet_eta = new THStack("jet_eta","");
        THStack* bjet_eta = new THStack("bjet_eta","");
        THStack* lepton_phi = new THStack("lepton_phi","");
        THStack* photon_phi = new THStack("photon_phi","");
        THStack* jet_phi = new THStack("jet_phi","");
        THStack* bjet_phi = new THStack("bjet_phi","");
        THStack* met = new THStack("met","");
	THStack* njets = new THStack("njets","");
	THStack* nbjets = new THStack("nbjets","");
	THStack* deltaR_p_l = new THStack("deltaR_p_l","");
	THStack* deltaR_p_j = new THStack("deltaR_p_j","");
	THStack* deltaR_p_b = new THStack("deltaR_p_b","");
	THStack* deltaR_l_j = new THStack("deltaR_l_j","");
        THStack* deltaR_l_b = new THStack("deltaR_l_b","");
        THStack* deltaR_j_b = new THStack("deltaR_j_b","");
	THStack* nleptons = new THStack("nleptons","");
	
	THStack* lepton_pt_reco = new THStack("lepton_pt_reco","");
        THStack* photon_pt_reco = new THStack("photon_pt_reco","");
        THStack* jet_pt_reco = new THStack("jet_pt_reco","");
        THStack* bjet_pt_reco = new THStack("bjet_pt_reco","");
        THStack* lepton_eta_reco = new THStack("lepton_eta_reco","");
        THStack* photon_eta_reco = new THStack("photon_eta_reco","");
        THStack* jet_eta_reco = new THStack("jet_eta_reco","");
        THStack* bjet_eta_reco = new THStack("bjet_eta_reco","");
        THStack* lepton_phi_reco = new THStack("lepton_phi_reco","");
        THStack* photon_phi_reco = new THStack("photon_phi_reco","");
        THStack* jet_phi_reco = new THStack("jet_phi_reco","");
        THStack* bjet_phi_reco = new THStack("bjet_phi_reco","");
        THStack* met_reco = new THStack("met_reco","");
	THStack* njets_reco = new THStack("njets_reco","");
	THStack* nbjets_reco = new THStack("nbjets_reco","");
	THStack* deltaR_p_l_reco = new THStack("deltaR_p_l_reco","");
	THStack* deltaR_p_j_reco = new THStack("deltaR_p_j_reco","");
	THStack* deltaR_p_b_reco = new THStack("deltaR_p_b_reco","");
	THStack* deltaR_l_j_reco = new THStack("deltaR_l_j_reco","");
        THStack* deltaR_l_b_reco = new THStack("deltaR_l_b_reco","");
        THStack* deltaR_j_b_reco = new THStack("deltaR_j_b_reco","");
	THStack* nleptons_reco = new THStack("nleptons_reco","");


	std::vector<THStack*> stacks = {bjet_eta, bjet_phi, bjet_pt, 
					deltaR_j_b, deltaR_l_b, deltaR_l_j, deltaR_p_b, deltaR_p_j, deltaR_p_l, 
					jet_eta, jet_phi, jet_pt, 
					lepton_eta, lepton_phi, lepton_pt, 
					photon_eta, photon_phi, photon_pt, 
					met, nbjets, njets,
					nleptons};
	
	std::vector<THStack*> stacks_reco = {bjet_eta_reco, bjet_phi_reco, bjet_pt_reco, 
					deltaR_j_b_reco, deltaR_l_b_reco, deltaR_l_j_reco, deltaR_p_b_reco, deltaR_p_j_reco, deltaR_p_l_reco, 
					jet_eta_reco, jet_phi_reco, jet_pt_reco, 
					lepton_eta_reco, lepton_phi_reco, lepton_pt_reco, 
					photon_eta_reco, photon_phi_reco, photon_pt_reco, 
					met_reco, nbjets_reco, njets_reco,
					nleptons_reco};



	std::vector<TString> png_names = {"bjet_eta.png", "bjet_phi.png", "bjet_pt.png", 
					"deltar_j_b.png", "deltar_l_b.png","deltar_l_j.png","deltar_p_b.png","deltar_p_j.png","deltar_p_l.png",
					"jet_eta.png", "jet_phi.png", "jet_pt.png", 
					"lepton_eta.png", "lepton_phi.png", "lepton_pt.png",
					"photon_eta.png", "photon_phi.png", "photon_pt.png",
					"met.png", "nbjets.png", "njets.png", "nleptons.png"};
	
	std::vector<TString> png_names_reco = {"bjet_eta_reco.png", "bjet_phi_reco.png", "bjet_pt_reco.png",
                                        "deltar_j_b_reco.png", "deltar_l_b_reco.png","deltar_l_j_reco.png","deltar_p_b_reco.png","deltar_p_j_reco.png","deltar_p_l_reco.png",
                                        "jet_eta_reco.png", "jet_phi_reco.png", "jet_pt_reco.png",
                                        "lepton_eta_reco.png", "lepton_phi_reco.png", "lepton_pt_reco.png",
                                        "photon_eta_reco.png", "photon_phi_reco.png", "photon_pt_reco.png",
                                        "met_reco.png", "nbjets_reco.png", "njets_reco.png","nleptons_reco.png"};

	
	std::vector<TString> hist_names = {"hBJetEta", "hBJetPhi", "hBJetPt",
					"deltaR_jb", "deltaR_lb", "deltaR_lj", "deltaR_pb", "deltaR_pj", "deltaR_pl",
				       	"hJetEta", "hJetPhi", "hJetPt",
					"hLeptonEta", "hLeptonPhi", "hLeptonPt",
					"hPhotonEta", "hPhotonPhi", "hPhotonPt",
					"hMET", "hNBJets", "hNJets", "hNLeptons"};	
	
	std::vector<int> bins = {50, 50, 150,
				50, 50, 50, 50, 50, 50, 50,
				50, 50, 150,
				50, 50, 150,
				50, 50, 150,
				150, 10, 10, 10};

	std::vector<TString> xaxis = {"Leading BJet Eta", "Leading BJet Phi (rad)", "Leading BJet pT (GeV)",
					"delta R between Jet and BJet", "delta R between Lepton and BJet", "delta R between Lepton and Jet", "delta R between Photon and BJet", "delta R between Jet and Photon" , "delta R between Photon and Lepton",
					"Leading Jet Eta", "Leading Jet Phi (rad)", "Leading Jet pT (GeV)",
					"Leading Lepton Eta", "Leading Lepton Phi (rad)", "Leading Lepton pT (GeV)",
					"Leading Photon Eta", "Leading Photon Phi (rad)", "Leading Photon pT (GeV)",
					"MET (GeV)", "Number of B-Jets", "Number of Jets", "Number of Leptons"};

	std::vector<int> max_x = {6, 6, 450, 
					10, 10, 10, 10, 10, 10, 
					6, 6, 450, 
					6, 6, 450,
					6, 6, 600,
					300, 10, 10, 10};

	std::vector<int> min_x;

	for (int i = 0; i < max_x.size(); i++){
		if( max_x[i] == 6 ){
			min_x.push_back(-6);
		}
		else{ min_x.push_back(0); }
	}

        //OPEN HIST FILES
        //TFile* f1 = new TFile("hadronic_distributions.hist");
        TFile* f2 = new TFile("dilepton_distributions_gen.hist");
        TFile* f3 = new TFile("singlelep_distributions_gen.hist");
        TFile* f4 = new TFile("signal_distributions_gen.hist");

	TFile* f5 = new TFile("dilepton_distributions_reco.hist");
        TFile* f6 = new TFile("singlelep_distributions_reco.hist");
        TFile* f7 = new TFile("signal_distributions_reco.hist");
	
	for( int i = 0; i < stacks_reco.size(); i++){

		delete gROOT->FindObject("h1");
		delete gROOT->FindObject("d1");
		delete gROOT->FindObject("s1");

		TH1F* h1 = new TH1F("h1", "Signal", bins[i], min_x[i] , max_x[i]);
	        h1 = (TH1F*)f7->Get(hist_names[i]);
        	h1->SetFillColor(kBlue);
        	h1->SetFillStyle(3001);
        	//lepton_pt->Add(h1, "nostack");

	        TH1F* d1 = new TH1F("d1", "Dilepton", bins[i], min_x[i] , max_x[i]);
        	d1 = (TH1F*)f5->Get(hist_names[i]);
        	d1->SetFillColor(kGreen);
        	stacks_reco[i]->Add(d1);

        	TH1F* s1 = new TH1F("s1", "Single Lepton", bins[i], min_x[i] , max_x[i]);
        	s1 = (TH1F*)f6->Get(hist_names[i]);
        	s1->SetFillColor(kYellow);
        	stacks_reco[i]->Add(s1);

        	TCanvas* x1 = new TCanvas("x1", "x1", 1000, 1000);
        	stacks_reco[i]->Draw("hist");
        	h1->SetDrawOption("hist");
        	h1->Draw("hist same");
        	stacks_reco[i]->GetXaxis()->SetTitle(xaxis[i]);
        	gPad->BuildLegend(0.75, 0.75, 0.95, 0.95, "");

        	x1->Print(png_names_reco[i]);

		if(x1) { x1->Close();gSystem->ProcessEvents(); }
	}
	
	for( int i = 0; i < stacks.size(); i++){

		delete gROOT->FindObject("h1");
		delete gROOT->FindObject("d1");
		delete gROOT->FindObject("s1");

		TH1F* h1 = new TH1F("h1", "Signal", bins[i], min_x[i] , max_x[i]);
	        h1 = (TH1F*)f4->Get(hist_names[i]);
        	h1->SetFillColor(kBlue);
        	h1->SetFillStyle(3001);
        	//lepton_pt->Add(h1, "nostack");

	        TH1F* d1 = new TH1F("d1", "Dilepton", bins[i], min_x[i] , max_x[i]);
        	d1 = (TH1F*)f2->Get(hist_names[i]);
        	d1->SetFillColor(kGreen);
        	stacks[i]->Add(d1);

        	TH1F* s1 = new TH1F("s1", "Single Lepton", bins[i], min_x[i] , max_x[i]);
        	s1 = (TH1F*)f3->Get(hist_names[i]);
        	s1->SetFillColor(kYellow);
        	stacks[i]->Add(s1);

        	TCanvas* x1 = new TCanvas("x1", "x1", 1000, 1000);
        	stacks[i]->Draw("hist");
        	h1->SetDrawOption("hist");
        	h1->Draw("hist same");
        	stacks[i]->GetXaxis()->SetTitle(xaxis[i]);
        	gPad->BuildLegend(0.75, 0.75, 0.95, 0.95, "");

        	x1->Print(png_names[i]);

		if(x1) { x1->Close();gSystem->ProcessEvents(); }
	}

		
	return 0;
	}
