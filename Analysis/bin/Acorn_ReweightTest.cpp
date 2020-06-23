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
#include "Acorn/NTupler/interface/Reduction.h"
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

  // StandaloneReweight rw_hel("Acorn.Analysis.standalone_reweight", "rw_WG-HEL");

  TFile fin(argc >= 2 ? argv[1] : "EventTree.root");
  std::string model = argc >= 3 ? argv[2] : "rw_WG-EWDim6";
  StandaloneReweight rw("Acorn.Analysis.standalone_reweight", model);
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
  for (unsigned i = 0; i < 30; ++i) {
  	evt.SetEvent(i);
    auto const& lheparts = evt.GetPtrVec<ac::GenParticle>("lheParticles");

  	std::vector<std::vector<double>> parts;
  	std::vector<int> pdgs;
  	std::vector<int> hels;
  	std::vector<int> stats;

    bool skip = false;

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
        skip = true;
        pdgs.push_back(replace_pdg[p->pdgId()]);
      } else {
    		pdgs.push_back(p->pdgId());
      }
  		hels.push_back(int(p->spin()));
  		stats.push_back(p->status());
  	}

    if (skip) continue;

    auto const* info = evt.GetPtr<ac::EventInfo>("eventInfo");
    bool verbose = false;
    auto wts = rw.ComputeWeights(parts, pdgs, hels, stats, info->lheAlphaS(), false, verbose);
    auto wts_trans = rw.TransformWeights(wts);
    int idx = 1;
    if (model == "rw_WG-HEL") {
      idx = 19;
    }
    std::cout << model << ": " << wts_trans[idx] << "\t" << wts_trans[idx + 1] << "\n";

    // auto wts_hel = rw_hel.ComputeWeights(parts, pdgs, hels, stats, info->lheAlphaS(), false, verbose);
    // auto wts_hel_trans = rw_hel.TransformWeights(wts_hel);
    // std::cout << "HEL: " << wts_hel_trans[19] << "\t" << wts_hel_trans[20] << "\n";

    // if (verbose) {
    //   for (unsigned i = 0; i < wts.size(); ++i) {
    //     std::cout << "[" << i << "] = " << wts[i] << "\n";
    //   }
    // }
    // double w1 = (info->lheWeights().at(100000) + 1.0);
    // double w2 = (info->lheWeights().at(100004) + 1.0);
    // auto wts_hel = rw.ComputeWeights(parts, pdgs, hels, stats, info->lheAlphaS(), true, verbose);

    // double emulate_w1 = wts[0] / wts[2];
    // double emulate_w2 = wts[1] / wts[2];

    // double round_w1 = reduceMantissaToNbitsRounding(emulate_w1 - 1.0, 10);
    // double round_w2 = reduceMantissaToNbitsRounding(emulate_w2 - 1.0, 10);
    // double round_nosub_w1 = reduceMantissaToNbitsRounding(emulate_w1, 10);
    // double round_nosub_w2 = reduceMantissaToNbitsRounding(emulate_w2, 10);
    // double no_round_w1 = emulate_w1 - 1.0;
    // double no_round_w2 = emulate_w2 - 1.0;

    // double relative = ( (wts[1] / wts[0]) / (w2 / w1)  );
    // if (std::abs(relative - 1.0) > 0.05) verbose = true;
    // if (verbose) {
    //   std::cout << "Emulate: " << emulate_w1 << "\t" << emulate_w2 << "\t" << (emulate_w2 / emulate_w1) << "\n";
    //   std::cout << "Round: " << round_w1 << "\t" << round_w2 << "\t" << ((round_w2 + 1.) / (round_w1 + 1.) ) << "\n";
    //   std::cout << "NoRound: " << no_round_w1 << "\t" << no_round_w2 << "\t" << ((no_round_w2 + 1.) / (no_round_w1 + 1.) ) << "\n";
    //   std::cout << "RoundNoSub: " << round_nosub_w1 << "\t" << round_nosub_w2 << "\t" << ((round_nosub_w2 ) / (round_nosub_w1) ) << "\n";
    //   std::cout << "Actual: " << info->lheWeights().at(100000) << "\t" << info->lheWeights().at(100004) << "\n";
    //   std::cout << (w2 / w1) << "\t" << (wts[1] / wts[0]) << "\t" << relative << "\n";
    // }

    // std::cout << "eventInfo[0]: " << (info->lheWeights().at(100000) + 1.0) << "\n";
    // std::cout << "eventInfo[1]: " << (info->lheWeights().at(100001) + 1.0) << "\n";
    // for (auto wt : wts) std::cout << wt << "\t";
    // std::cout << "\n";
  }
}
