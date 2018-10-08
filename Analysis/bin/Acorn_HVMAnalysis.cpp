#include <iostream>
#include <vector>
#include <string>
#include <fstream>
// Services
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/AnalysisBase.h"
#include "Acorn/Analysis/interface/AnalysisSequence.h"
#include "Acorn/NTupler/interface/json.hpp"
#include "Acorn/NTupler/interface/Reduction.h"
// Modules
#include "Acorn/Analysis/interface/GenericModule.h"
#include "Acorn/Analysis/interface/HVMGenAnalysis.h"
#include "Acorn/Analysis/interface/EventCounters.h"
#include "Acorn/Analysis/interface/LumiMask.h"
#include "Acorn/Analysis/interface/SampleStitching.h"
#include "Compression.h"
using std::string;
using std::vector;
using std::set;

std::vector<std::string> ParseFileLines(std::string const& file_name) {
  // Build a vector of input files
  std::vector<std::string> files;
  std::ifstream file;
  file.open(file_name.c_str());
  if (!file.is_open()) {
    std::cerr << "Warning: File " << file_name << " cannot be opened." << std::endl;
    return files;
  }
  std::string line = "";
  while (std::getline(file, line)) {  // while loop through lines
    files.push_back(line);
  }
  file.close();
  return files;
}

std::vector<std::string> GetFilesForJob(std::vector<std::string> const& filelists, std::string const& prefix, unsigned offset, unsigned step) {
  vector<string> files;
  for (auto const& filelist : filelists) {
    auto i_files = ParseFileLines(filelist);
    files.insert(files.end(), i_files.begin(), i_files.end());
  }
  vector<string> do_files;
  for (auto & f : files) {
    f = prefix + f;
  }

  for (unsigned i = offset; i < files.size(); i += step) {
    do_files.push_back(files[i]);
  }
  return do_files;
}

template<typename T, typename Range>
bool contains(Range const& r, T const& value)
{
  return std::find(r.begin(), r.end(), value) != r.end();
}



int main(int argc, char* argv[]) {
  using json = nlohmann::json;

  json js = json::parse(argv[1]);
  json const& jsc = js;

  // Should move to passing the list of files to process directly?
  //std::string outname = "testout.root";
  std::string outname = jsc["output"];
  vector<string> do_files = GetFilesForJob(jsc["filelists"], "", jsc["file_offset"], jsc["file_step"]);
  //vector<string> do_files;
  //do_files.push_back("file:///nfs/dust/cms/user/dewita/CMSSW_9_4_4/src/Acorn/NTupler/test/EventTree_omega.root");
  //std::string outputdir = "output/";

  std::string outputdir = jsc["outdir"];
  

  std::set<std::string> sequences = jsc["sequences"];
  std::map<std::string, std::shared_ptr<fwlite::TFileService>> fs;
 // fs["HVMGen"] = std::make_shared<fwlite::TFileService>(
  //    outputdir + "/" + "HVMGen" + "/" + outname);
  //fs["HVMGen"]->file().SetCompressionSettings(ROOT::CompressionSettings(ROOT::kZLIB, 5));
  for (auto const& seq : sequences) {
    system(("mkdir -p "+outputdir+"/"+seq).c_str());
    fs[seq] = std::make_shared<fwlite::TFileService>(
        outputdir + "/" + seq + "/" + outname);
    fs[seq]->file().SetCompressionSettings(ROOT::CompressionSettings(ROOT::kZLIB, 5));
  }


  //ac::AnalysisBase analysis("HVMAnalysis", do_files, "EventTree", 1000);
  ac::AnalysisBase analysis("HVMAnalysis", do_files, "EventTree", jsc.count("events") ? int(jsc["events"]) : -1);
  analysis.SetTTreeCaching(true);
  analysis.StopOnFileFailure(true);
  analysis.RetryFileAfterFailure(7, 3);
  analysis.CalculateTimings(false);

  //bool is_data = contains(jsc["attributes"], "data");

  ac::Sequence hvmgen_seq;
  if (sequences.count("HVMGen")) {
    auto hvmgen_fs = fs.at("HVMGen").get();
    hvmgen_seq.BuildModule(ac::EventCounters("EventCounters").set_fs(hvmgen_fs));
 
    hvmgen_seq.BuildModule(ac::HVMGenAnalysis("HVMGenAnalysis")
                             .set_fs(hvmgen_fs));

    hvmgen_seq.InsertSequence("HVMGen", analysis);
  }

  analysis.RunAnalysis();
  return 0;
}

