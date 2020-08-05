#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <vector>
#include <string>
#include <iostream>

#include "Acorn/NTupler/interface/GenParticle.h"

int BasicAnalysis() {
	std::vector<std::string> inputFiles = {
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_1.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_10.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_11.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_12.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_13.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_14.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_15.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_16.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_17.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_2.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_3.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_4.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_5.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_6.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_7.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_8.root",
		"root://eoscms.cern.ch//store/group/phys_smp/agilbert/wgamma_2018_v5/TGJets_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/TGJets_leptonDecays-amcatnlo/200427_185831/0000/EventTree_9.root"
	};

	unsigned maxEvents = 10;
	unsigned iEvt = 0;
	for (auto const& fileName : inputFiles) {
		TFile *f = TFile::Open(fileName.c_str());


		TTreeReader reader("EventTree", f);

		TTreeReaderValue<std::vector<ac::GenParticle>> genParts(reader, "genParticles");
		while (reader.Next()) {
			std::cout << "--- Event ---\n";
			ac::GenParticle const* lepton = nullptr;
			ac::GenParticle const* neutrino = nullptr;
			for (auto const& part : *genParts) {
				part.Print();
				if (part.status() == 1 && (abs(part.pdgId()) == 11 || abs(part.pdgId()) == 13)) {
					// If we already found a lepton, only take this one if it has higher pT
					if (!lepton || part.pt() > lepton->pt()) {
						lepton = &part;
					}
				}
				if (part.status() == 1 && (abs(part.pdgId()) == 12 || abs(part.pdgId()) == 14)) {
					// If we already found a neutrino, only take this one if it has higher pT
					if (!neutrino || part.pt() > neutrino->pt()) {
						neutrino = &part;
					}
				}
			}
			++iEvt;

			// Try to find a W boson...
			if (lepton && neutrino) {
				std::cout << "--> W boson with mass " << (lepton->vector() + neutrino->vector()).M() << "\n";
			}

			if (iEvt >= maxEvents) break;
		}
		f->Close();
		if (iEvt >= maxEvents) break;
	}
	return 0;
}