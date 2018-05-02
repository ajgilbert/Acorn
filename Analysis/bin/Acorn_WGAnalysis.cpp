#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/AnalysisBase.h"
#include "Acorn/Analysis/interface/AnalysisSequence.h"
#include "Acorn/Analysis/interface/GenericModule.h"
#include "Acorn/Analysis/interface/WGAnalysis.h"
using std::string;
using std::vector;
using std::set;

int main(int argc, char* argv[]) {
  vector<string> do_files = {"EventTree.root"};

  std::string outputdir = "./output";
  std::map<std::string, std::shared_ptr<fwlite::TFileService>> fs;
  for (auto const& seq : {"Main"}) {
    fs[seq] = std::make_shared<fwlite::TFileService>(
        outputdir + "/" + seq + "/" + "output.root");
  }

  ac::AnalysisBase analysis("WGAnalysis", do_files, "EventTree", -1);
  analysis.SetTTreeCaching(true);
  analysis.StopOnFileFailure(true);
  analysis.RetryFileAfterFailure(7, 3);
  analysis.CalculateTimings(true);


  ac::Sequence main_seq;
  main_seq.BuildModule(ac::WGAnalysis("WGAnalysis").set_fs(fs.at("Main").get()));

  main_seq.InsertSequence("Main", analysis);

  analysis.RunAnalysis();
}
