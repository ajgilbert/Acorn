#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <vector>
#include <string>
#include <iostream>


#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/GenParticle.h"

namespace { //fields for interactive configuration
  int verbose_ = 0; //how verbose the printout info should be
  unsigned maxEvents_ = 0; //0 means process all events
  bool fillAllEvents_ = true; //fill histograms for basic selection or all events
  TString outfile_ = "BasicAnalysisTwo.hist"; //name of output file
}

// check if a jet is a b-jet by geometric matching
bool isBJet(ac::Candidate const& jet, ac::GenParticle const& bquark) {
  TLorentzVector lv1, lv2;
  lv1.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.M());
  lv2.SetPtEtaPhiM(bquark.pt(), bquark.eta(), bquark.phi(), bquark.M());
  double deltaRMax = 0.3;
  bool retval = abs(lv1.DeltaR(lv2)) < deltaRMax;
  if(verbose_ > 2) std::cout << "isBJet check = " << retval << std::endl;
  return retval;
}

int BasicAnalysisTwo() {
  //list of files to process
  std::vector<std::string> inputFiles = { "/home/mmackenz/storage/TGJets/EventTree_1.root"};

  //output file
  TFile* fOut = new TFile(outfile_.Data(), "RECREATE");
  fOut->cd();
  
  //histograms to fill
  TH1F* hPhotonPt = new TH1F("hPhotonPt", "Leading Photon pT"      , 100, 0., 200.);
  TH1F* hLeptonPt = new TH1F("hLeptonPt", "Leading Lepton pT"      , 100, 0., 200.);
  TH1F* hBJetPt   = new TH1F("hBJetPt"  , "Leading b-tagged Jet pT", 100, 0., 200.);
  TH1F* hJetPt    = new TH1F("hJetPt"   , "Leading Jet pT"         , 100, 0., 200.);
  TH1F* hNBJets   = new TH1F("hNBJets"  , "N(b-tagged Jets)"       ,  10, 0.,  10.);
  TH1F* hNJets    = new TH1F("hNJets"   , "N(untagged Jets)"       ,  10, 0.,  10.);
  TH1F* hMET      = new TH1F("hMET"     , "MET"                    , 100, 0., 200.);
  
  unsigned iEvt = 0; //track event number
  for (auto const& fileName : inputFiles) {
    TFile *f = TFile::Open(fileName.c_str(), "READ");
    TTreeReader reader("EventTree", f);
    TTreeReaderValue<std::vector<ac::GenParticle>> genParts(reader, "genParticles");
    TTreeReaderValue<std::vector<ac::Candidate>> genJets(reader, "genJets");

    while (reader.Next()) {
      if(iEvt % 10000 == 0) std::cout << "Processing event " << iEvt+1 << "...\n";
      ++iEvt;      
      if(verbose_ > 0) std::cout << "--- Event " << iEvt << "---\n";
      ac::GenParticle const* lepton = nullptr;
      ac::GenParticle const* neutrino = nullptr;
      ROOT::Math::PtEtaPhiMVector met;
      std::vector<ac::GenParticle> bquarks;
      std::vector<ac::GenParticle> photons;
      for (auto const& part : *genParts) {
	if(verbose_ > 1) part.Print();
	//look for muons and electrons
	if (part.status() == 1 && (abs(part.pdgId()) == 11 || abs(part.pdgId()) == 13)) {
	  // If we already found a lepton, only take this one if it has higher pT
	  if (!lepton || part.pt() > lepton->pt()) {
	    lepton = &part;
	  }
	}
	//look for neutrinos
	if (part.status() == 1 && (abs(part.pdgId()) == 12 || abs(part.pdgId()) == 14 || abs(part.pdgId()) == 16)) {
	  // If we already found a neutrino, only take this one if it has higher pT
	  if (!neutrino || part.pt() > neutrino->pt()) {
	    neutrino = &part;
	  }
	  //add the neutrino to the MET
	  met += neutrino->vector();
	}
	//look for b-quarks
	if(abs(part.pdgId()) == 5) {
	  bquarks.push_back(part);
	}
	//look for photons
	if(part.pdgId() == 22) {
	  photons.push_back(part);
	}
      } //end gen particle loop
      
      //loop through jets and identify b-jets
      std::vector<ac::Candidate> bjets;
      std::vector<ac::Candidate> jets;
      for (ac::Candidate const& jet : *genJets) {
	if(verbose_ > 1) {std::cout << "jet: "; jet.Print();}
	bool isBTagged = false;
	for(auto const& bquark : bquarks) {
	  if(isBJet(jet, bquark)) {
	    isBTagged = true;
	    break;
	  }
	}	
	if(isBTagged)
	  bjets.push_back(jet);
	else
	  jets.push_back(jet);	
      }
      
      //highest pT objects
      ac::Candidate const* bjet = nullptr;
      ac::Candidate const* jet = nullptr;
      ac::GenParticle const* photon = nullptr;
      //look for highest pT photon
      for(ac::GenParticle const& p : photons) {
      	if(!photon || photon->pt() < p.pt())
      	  photon = &p;
      }
      //look for highest pT bjet
      for(ac::Candidate const& b : bjets) {
      	if(!bjet || bjet->pt() < b.pt())
      	  bjet = &b;
      }
      //look for highest pT untagged jet
      for(ac::Candidate const& j : jets) {
      	if(!jet || jet->pt() < j.pt())
      	  jet = &j;
      }
      //basic event selection, require each object is found
      if(bjet && jet && photon && lepton) {
	if(verbose_ > 0) {
	  std::cout << "Accepted event!:\n";	
	  std::cout << "lepton = "; photon->Print();
	  std::cout << "photon = "; lepton->Print();
	  std::cout << "bjet   = "; bjet->Print();
	  std::cout << "jet    = "; jet->Print();
	}
	//fill histograms
	hPhotonPt->Fill(photon->pt());
	hLeptonPt->Fill(lepton->pt());
	hBJetPt  ->Fill(bjet->pt());
	hJetPt   ->Fill(jet->pt());
	hNBJets  ->Fill(bjets.size());
	hNJets   ->Fill(jets.size());
	hMET     ->Fill(met.Pt());
      } else if(fillAllEvents_) {
	hPhotonPt->Fill(((photon) ? photon->pt() : -1.));
	hLeptonPt->Fill(((lepton) ? lepton->pt() : -1.));
	hBJetPt  ->Fill(((bjet) ? bjet->pt() : -1.));
	hJetPt   ->Fill(((jet) ? jet->pt() : -1.));
	hNBJets  ->Fill(bjets.size());
	hNJets   ->Fill(jets.size());
	hMET     ->Fill(met.Pt());
      }
	
      // Try to find a W boson...
      if (lepton && neutrino) {
	if(verbose_ > 0) std::cout << "--> W boson with mass " << (lepton->vector() + neutrino->vector()).M() << "\n";
      }
      //print object counts if verbose output
      if(verbose_ > 0) std::cout << "N(leptons) = " << (!lepton ? 0 : 1) << std::endl
				 << "N(photons) = " << photons.size() << std::endl
				 << "N(bjets  ) = " << bjets.size() << std::endl
				 << "N(jets   ) = " << jets.size() << std::endl;
      if (maxEvents_ > 0 && iEvt >= maxEvents_) break;
    }
    f->Close();
    if (maxEvents_ > 0 && iEvt >= maxEvents_) break;
  }
  //write the histograms to disk
  fOut->cd();
  hPhotonPt->Write();
  hLeptonPt->Write();
  hBJetPt->Write();
  hJetPt->Write();
  hNBJets->Write();
  hNJets->Write();
  fOut->Write();
  fOut->Close();
  
  return 0;
}
