#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
// Services
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/AnalysisBase.h"
#include "Acorn/Analysis/interface/AnalysisSequence.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/json.hpp"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/GenParticle.h"
// Modules
#include "Acorn/Analysis/interface/GenericModule.h"
#include "Acorn/Analysis/interface/WGAnalysis.h"
#include "Acorn/Analysis/interface/WGDataAnalysis.h"
#include "Acorn/Analysis/interface/JESStudy.h"
#include "Acorn/Analysis/interface/WGAnalysisTools.h"
#include "Acorn/Analysis/interface/WGTagAndProbe.h"
#include "Acorn/Analysis/interface/DiMuonAnalysis.h"
#include "Acorn/Analysis/interface/EventCounters.h"
#include "Acorn/Analysis/interface/LumiMask.h"
#include "Acorn/Analysis/interface/SampleStitching.h"
#include "Acorn/Analysis/interface/StandaloneReweight.h"
#include "Compression.h"
using std::string;
using std::vector;
using std::set;


int main(int argc, char* argv[]) {

  StandaloneReweight rw("Acorn.Analysis.standalone_reweight", "rw_WG-EWDim6");

  TFile fin("EventTree.root");
  TTree *tin = (TTree*)fin.Get("EventTree");
  ac::TreeEvent evt;
  evt.SetTree(tin);
  std::map<int, int> replace_pdg = {
    {11, 13},
    {12, 14},
    {15, 13},
    {16, 14},
    {-11, -13},
    {-12, -14},
    {-15, -13},
    {-16, -14}
  };
  for (unsigned i = 0; i < tin->GetEntries(); ++i) {
  	evt.SetEvent(i);
    auto const& lheparts = evt.GetPtrVec<ac::GenParticle>("lheParticles");

  	std::vector<std::vector<double>> parts;
  	std::vector<int> pdgs;
  	std::vector<int> hels;
  	std::vector<int> stats;

  	for (auto const& p : lheparts) {
  		if (std::abs(p->status()) != 1) {
  			continue;
  		}
  		if (p->status() == -1) {
	  		parts.push_back({p->M(), 0., 0., p->pt()});
  		}
  		if (p->status() == 1) {
  			parts.push_back({p->energy(), p->vector().px(), p->vector().py(), p->vector().pz()});
  		}
      if (replace_pdg.count(p->pdgId())) {
        pdgs.push_back(replace_pdg[p->pdgId()]);
      } else {
    		pdgs.push_back(p->pdgId());
      }
  		hels.push_back(int(p->spin()));
  		stats.push_back(p->status());
  	}

    auto const* info = evt.GetPtr<ac::EventInfo>("eventInfo");

    auto wts = rw.ComputeWeights(parts, pdgs, hels, stats, info->lheAlphaS(), false, true);
    double w1 = (info->lheWeights().at(100000) + 1.0);
    double w2 = (info->lheWeights().at(100004) + 1.0);
    auto wts_hel = rw.ComputeWeights(parts, pdgs, hels, stats, info->lheAlphaS(), true, true);

    std::cout << (w2 / w1) << "\t" << (wts[1] / wts[0]) << "\t" << (wts_hel[1] / wts_hel[0]) << "\n";
    // std::cout << "eventInfo[0]: " << (info->lheWeights().at(100000) + 1.0) << "\n";
    // std::cout << "eventInfo[1]: " << (info->lheWeights().at(100001) + 1.0) << "\n";
    // for (auto wt : wts) std::cout << wt << "\t";
    // std::cout << "\n";
  }
}
