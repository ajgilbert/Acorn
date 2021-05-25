void jet_analysis(){

		TH1F* NJets1 = new TH1F("NJet1", "Signal", 10, 0, 10);
	        TH1F* NJets2 = new TH1F("NJet2", "Signal", 10, 0, 10);

		int total_jets = 0;
		int cuts_lev = 0;
		int photon_count =0;
		int lepton_count = 0;
		int bjet_count = 0;
		int iEvt = 0;
		std::vector<std::string> input_files =   {"/home/mmackenz/storage/TGJets/EventTree_1.root",//};//,
                                                        "/home/mmackenz/storage/TGJets/EventTree_2.root",
                                                        //"/home/mmackenz/storage/TGJets/EventTree_3.root",
                                                        //"/home/mmackenz/storage/TGJets/EventTree_4.root",
                                                        "/home/mmackenz/storage/TGJets/EventTree_5.root"};
		double event_count;
		double pT_gamma = 50;
	        double eta_gamma = 2.5;
        	double ET_miss = 20;
        	double pT_lepton = 20;
        	double eta_lepton = 2.5;
        	double pT_bjet = 20;
        	double eta_b = 2.5;
        	double eta_j = 5.0;

	
			
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
				
				int jet_check = 0;
				int photon_check = 0;
				int lepton_check = 0;
				int bjet_check = 0;
				int jetsize = jets.size();

				if(jet && met->pt()>ET_miss){
					
						total_jets += eventInfo->totalWeight();
						for(auto const& j : jets){
							//j.Print();
							for (auto const& k : jets){
								double overlap = sqrt(pow((j.eta() - k.eta()),2) + pow((j.phi() - k.phi()),2));
								if(overlap != 0 && overlap < 0.5){
									jet_check = 1;
									//j.Print();
									//k.Print();
								}
							}
							for (auto const& k : photons){
                                                                double overlap = sqrt(pow((j.eta() - k.eta()),2) + pow((j.phi() - k.phi()),2));
                                                                if(overlap != 0 && overlap < 0.5){
                                                                        photon_check = 1;
                                                                        //j.Print();
                                                                        //k.Print();
                                                                }
                                                        }
							for (auto const& k : leptons){
                                                                double overlap = sqrt(pow((j.eta() - k.eta()),2) + pow((j.phi() - k.phi()),2));
                                                                if(overlap != 0 && overlap < 0.5){
                                                                        lepton_check = 1;
                                                                        //j.Print();
                                                                        //k.Print();
                                                                }
                                                        }
							for (auto const& k : bjets){
                                                                double overlap = sqrt(pow((j.eta() - k.eta()),2) + pow((j.phi() - k.phi()),2));
                                                                if(overlap != 0 && overlap < 0.5){
                                                                        bjet_check = 1;
                                                                        //j.Print();
                                                                        //k.Print();
                                                                }
                                                        }



						}
						if(photon_check == 1){
							photon_count += eventInfo->totalWeight();
						}
						if(lepton_check == 1){
                                                        lepton_count += eventInfo->totalWeight();
                                                }
						if(bjet_check == 1){
                                                        bjet_count += eventInfo->totalWeight();
                                                }


						if(jet_check ==1){
							jetsize = jets.size() - 1;
							event_count += eventInfo->totalWeight();
							for(auto const& j : jets){
                                                        	j.Print();
							}		
						cout<<"----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
					
						}
					
				NJets1->Fill(jets.size(), eventInfo->totalWeight());
				NJets2->Fill(jetsize, eventInfo->totalWeight());
				}

		}
}
cout<<"Overlapping Jets = "<< event_count/2<<"\n"<<endl;
cout<<"Overlapping Photon = "<< photon_count<<"\n"<<endl;
cout<<"Overlapping Lepton = "<<lepton_count<< "\n"<<endl;
cout<<"Overlapping BJet = " <<bjet_count<<"\n"<<endl;
cout<<"Total Events = "<< total_jets<< "\n"<<endl;
TCanvas* x1 = new TCanvas("x1", "x1", 1000, 1000);
NJets1->SetFillColor(kGreen);
NJets1->Draw("hist");

TCanvas* x2 = new TCanvas("x2", "x2", 1000, 1000);
NJets2->SetFillColor(kBlue);
NJets2->Draw("hist");
}
