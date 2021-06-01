#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <vector>
#include <string>
#include <iostream>

#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/StandaloneReweight.h"

int EFTAnalysis() {

  StandaloneReweight rw("standalone_reweight", "rw_tG-SMEFTsim");
	std::vector<std::string> inputFiles = {
					       "/eos/cms/store/group/phys_smp/agilbert/wgamma_2018_test/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200819_145921/0000/EventTree_1.root"
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_1.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_10.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_11.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_12.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_13.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_14.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_15.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_16.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_17.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_2.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_3.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_4.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_5.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_6.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_7.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_8.root",
		// "root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_9.root"
	};

	unsigned maxEvents = 10;
	unsigned iEvt = 0;
	for (auto const& fileName : inputFiles) {
		TFile *f = TFile::Open(fileName.c_str());


		TTreeReader reader("EventTree", f);

		TTreeReaderValue<std::vector<ac::GenParticle>> genParts(reader, "lheParticles");
		TTreeReaderValue<ac::EventInfo> info(reader, "eventInfo");
		while (reader.Next()) {
			std::cout << "--- Event ---\n";

      std::vector<std::vector<double>> rw_parts;
      std::vector<int> pdgs;
      std::vector<int> hels;
      std::vector<int> stats;

      int top_idx = -1;
      for (unsigned ip = 0; ip < genParts->size(); ++ip) {
        bool keep = true;
        auto const& p = (*genParts)[ip];
        p.Print();
        int st = p.status();
        
        if (top_idx >= 0 && (ip - top_idx) <= 4) keep = false; // skip top decay products
        
        if (p.status() == 2 && std::abs(p.pdgId()) == 6) {
          top_idx = ip;
          st = 1; // treat the top as outgoing
        }
        
        if (!keep) continue; 

        if (st == -1) {
          rw_parts.push_back({p.M(), 0., 0., p.pt()});
        }
        if (st == 1) {
          rw_parts.push_back({p.energy(), p.vector().px(), p.vector().py(), p.vector().pz()});
        }
        pdgs.push_back(p.pdgId());
        hels.push_back(int(p.spin()));
        stats.push_back(st);
      }

      auto wts = rw.ComputeWeights(rw_parts, pdgs, hels, stats, info->lheAlphaS(), false, false);
      // std::cout << wts.size() << "\n";
      auto trans_wts = rw.TransformWeights(wts);
      for (auto wt : trans_wts) {
        std::cout << " - " << wt << "\n";
      }
      
			++iEvt;


			if (iEvt >= maxEvents) break;
		}
		f->Close();
		if (iEvt >= maxEvents) break;
	}
	return 0;
}
