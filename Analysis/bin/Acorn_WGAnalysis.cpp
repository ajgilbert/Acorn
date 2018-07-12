#include <iostream>
#include <vector>
#include <string>
#include <fstream>
// Services
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/AnalysisBase.h"
#include "Acorn/Analysis/interface/AnalysisSequence.h"
#include "Acorn/NTupler/interface/json.hpp"
// Modules
#include "Acorn/Analysis/interface/GenericModule.h"
#include "Acorn/Analysis/interface/WGAnalysis.h"
#include "Acorn/Analysis/interface/DiMuonAnalysis.h"
#include "Acorn/Analysis/interface/EventCounters.h"
#include "Acorn/Analysis/interface/LumiMask.h"
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
  std::string outname = jsc["output"];
  vector<string> do_files = GetFilesForJob(jsc["filelists"], "", jsc["file_offset"], jsc["file_step"]);
  std::string outputdir = jsc["outdir"];

  std::set<std::string> sequences = jsc["sequences"];
  std::map<std::string, std::shared_ptr<fwlite::TFileService>> fs;
  for (auto const& seq : sequences) {
    fs[seq] = std::make_shared<fwlite::TFileService>(
        outputdir + "/" + seq + "/" + outname);
  }

  ac::AnalysisBase analysis("DiMuonAnalysis", do_files, "EventTree", -1);
  analysis.SetTTreeCaching(true);
  analysis.StopOnFileFailure(true);
  analysis.RetryFileAfterFailure(7, 3);
  analysis.CalculateTimings(false);

  bool is_data = contains(jsc["attributes"], "data");

  ac::Sequence main_seq;
  if (sequences.count("Main")) {
    main_seq.BuildModule(ac::EventCounters("EventCounters").set_fs(fs.at("Main").get()));
    if (is_data) {
      main_seq.BuildModule(
          ac::LumiMask("LumiMask").set_fs(fs.at("Main").get()).set_input_file(jsc["data_json"]));
    }
    main_seq.BuildModule(ac::DiMuonAnalysis("DiMuonAnalysis")
                             .set_fs(fs.at("Main").get())
                             .set_year(jsc["year"])
                             .set_corrections("input/wgamma_corrections_2016_v1.root")
                             .set_is_data(is_data));

    main_seq.InsertSequence("Main", analysis);
  }

  ac::Sequence wg_gen_seq;

  if (sequences.count("wg_gen")) {
    wg_gen_seq.BuildModule(ac::WGAnalysis("WGAnalysis").set_fs(fs.at("wg_gen").get()));
    wg_gen_seq.InsertSequence("wg_gen", analysis);
  }

  analysis.RunAnalysis();
}
