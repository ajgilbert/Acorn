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


int main(int argc, char* argv[]) {
  using json = nlohmann::json;

  json js = json::parse(argv[1]);
  json const& jsc = js;
  std::string s_year = boost::lexical_cast<std::string>(jsc["year"]);
  int year = jsc["year"];

  // Should move to passing the list of files to process directly?
  // std::string outname = jsc["output"];
  vector<string> do_files = GetFilesForJob(jsc["filelists"], "", jsc["file_offset"], jsc["file_step"]);
  std::string outputdir = jsc["outdir"];

  std::set<std::string> sequences = jsc["sequences"];
  std::map<std::string, std::shared_ptr<fwlite::TFileService>> fs;
  for (auto const& seq : sequences) {
    std::string outfile = jsc["output_" + seq];
    fs[seq] = std::make_shared<fwlite::TFileService>(
        outputdir + "/" + outfile);
    fs[seq]->file().SetCompressionSettings(ROOT::CompressionSettings(ROOT::kZLIB, 5));
  }

  ac::AnalysisBase analysis("WGammaAnalysis", do_files, "EventTree", jsc.count("events") ? int(jsc["events"]) : -1);
  analysis.SetTTreeCaching(true);
  analysis.StopOnFileFailure(true);
  analysis.RetryFileAfterFailure(7, 3);
  analysis.CalculateTimings(false);

  bool is_data = ac::contains(jsc["attributes"], "data");

  ac::Sequence lumi_seq;
  if (sequences.count("Lumi")) {
    lumi_seq.BuildModule(ac::EventCounters("EventCounters").set_fs(fs.at("Lumi").get()));
    if (is_data) {
      lumi_seq.BuildModule(
          ac::LumiMask("LumiMask").set_fs(fs.at("Lumi").get()).set_input_file(jsc["data_json"]));
    }
    lumi_seq.InsertSequence("Lumi", analysis);
  }

  ac::Sequence dimuon_seq;
  if (sequences.count("DiMuon")) {
    dimuon_seq.BuildModule(ac::EventCounters("EventCounters").set_fs(fs.at("DiMuon").get()));
    if (is_data) {
      dimuon_seq.BuildModule(
          ac::LumiMask("LumiMask").set_fs(fs.at("DiMuon").get()).set_input_file(jsc["data_json"]));
    }
    dimuon_seq.BuildModule(ac::DiMuonAnalysis("DiMuonAnalysis")
                             .set_fs(fs.at("DiMuon").get())
                             .set_year(jsc["year"])
                             .set_corrections("wgamma/inputs/wgamma_corrections_2016_v3.root")
                             .set_is_data(is_data));

    dimuon_seq.InsertSequence("DiMuon", analysis);
  }

  std::map<std::string, ac::Sequence> wgamma_seqs;

  for (auto const& seq : sequences) {
    std::vector<std::string> as_vec;
    boost::split(as_vec, seq, boost::is_any_of("_"));
    if (!(as_vec.size() >= 1 && as_vec[0] == "WGamma")) continue;
    std::string subseq;
    if (as_vec.size() >= 2) subseq = as_vec[1];

    ac::Sequence & wgamma_seq = wgamma_seqs[seq];
    std::string wgamma_label = seq;
    auto wgamma_fs = fs.at(wgamma_label).get();

    std::vector<std::string> userDoubleNames = {"fixedGridRhoFastjetAll"};
    if (!is_data && (year == 2016 || year == 2017)) {
        userDoubleNames.push_back("NonPrefiringProb");
        userDoubleNames.push_back("NonPrefiringProbUp");
        userDoubleNames.push_back("NonPrefiringProbDown");
    }
    wgamma_seq.BuildModule(ac::GenericModule("UserDoubleUnpacker").set_function([=](ac::TreeEvent* event) {
        auto const& vals = event->GetPtr<ac::EventInfo>("eventInfo")->userDoubles();
        for (unsigned i = 0; i < userDoubleNames.size(); ++i) {
            event->Add(userDoubleNames[i], vals.at(i));
        }
        return 0;
    }));

    if (jsc.count("stitching")) {
      wgamma_seq.BuildModule(ac::SampleStitching("SampleStitching", jsc["stitching"]));
      // wgamma_seq.BuildModule(ac::SampleStitching("SampleStitching", jsc["stitching"]).set_fs(wgamma_fs));
    }


    int scale_weights = ac::ReadAttrValue<int>(jsc["attributes"], "scale_weights");
    if (scale_weights > 0) {
      wgamma_seq.BuildModule(ac::GenericModule("ScaleWeightUnpacker").set_function([=](ac::TreeEvent* event) {
          auto const* info = event->GetPtr<ac::EventInfo>("eventInfo");
          event->Add("scale_weights", ac::ExtractScaleVariations(*info, scale_weights));
          return 0;
      }));
    }

    int pdf_weights = ac::ReadAttrValue<int>(jsc["attributes"], "pdf_weights");
    std::vector<int> pdf_begin;
    std::vector<int> pdf_end;
    if (pdf_weights > 0 && subseq == "") { // only if the main sequence
      if (pdf_weights == 1) {
        // NNPDF31_nnl0_as_0118 + 100 replicas
        // NNPDF30_nlo_nf_5_pdfas + 100 replicas
        pdf_begin = {1114, 1002};
        pdf_end = {1214, 1102};
      }
    }

    int ps_weights = ac::ReadAttrValue<int>(jsc["attributes"], "ps_weights");

    auto counters = ac::EventCounters("EventCounters").set_fs(wgamma_fs);
    if (scale_weights > 0) {
      counters.AddWeightSet("scale_weights", 6, true);
    }
    wgamma_seq.BuildModule(counters);

    if (is_data) {
      wgamma_seq.BuildModule(
          ac::LumiMask("LumiMask").set_fs(wgamma_fs).set_input_file(jsc["data_json"]));
    }

    // Keep baseline set of vars
    int var_set = 1;
    // Unless this is some systematic variation - in which case keep only the core ones
    bool only_wg = true;
    if (subseq != "") {
      var_set = 0;
      only_wg = true;
    }

    int correct_p_energy = 0;
    int correct_e_energy = 0;
    int correct_m_energy = 0;
    int shift_met = 0;
    if (subseq == "EGMNoCorr") {
      correct_e_energy = -1;
      correct_p_energy = -1;
      correct_m_energy = -1;
    }
    if (subseq == "EScaleHi") correct_e_energy = 1;
    if (subseq == "EScaleLo") correct_e_energy = 2;
    if (subseq == "PScaleHi") correct_p_energy = 1;
    if (subseq == "PScaleLo") correct_p_energy = 2;
    if (subseq == "MScaleHi") correct_m_energy = 1;
    if (subseq == "MScaleLo") correct_m_energy = 2;
    if (subseq == "MetJesHi") shift_met = 1;
    if (subseq == "MetJesLo") shift_met = 2;
    if (subseq == "MetUncHi") shift_met = 3;
    if (subseq == "MetUncLo") shift_met = 4;

    std::cout << ac::ReadAttrValue<int>(jsc["attributes"], "scale_weights") << "\n";

    wgamma_seq.BuildModule(ac::WGDataAnalysis("WGDataAnalysis")
                             .set_fs(wgamma_fs)
                             .set_year(jsc["year"])
                             .set_corrections("wgamma/inputs/wgamma_corrections_" + s_year + "_v12.root")
                             .set_is_data(is_data)
                             .set_gen_classify("")
                             .set_do_wg_gen_vars(ac::contains(jsc["attributes"], "do_wg_gen_vars"))
                             .set_check_is_zg(ac::contains(jsc["attributes"], "check_is_zg"))
                             .set_check_is_wwg(ac::contains(jsc["attributes"], "check_is_wwg"))
                             .set_check_gen_mll(ac::contains(jsc["attributes"], "check_gen_mll"))
                             .set_do_presel(!ac::contains(jsc["attributes"], "no_presel"))
                             .set_only_wg(only_wg)
                             .set_correct_e_energy(correct_e_energy)
                             .set_correct_p_energy(correct_p_energy)
                             .set_correct_m_energy(correct_m_energy)
                             .set_shift_met(shift_met)
                             .set_rc_file("wgamma/inputs/muons/RoccoR" + s_year + ".txt")
                             .set_scale_weights(scale_weights)
                             .set_ps_weights(ps_weights)
                             .set_pdf_begin(pdf_begin)
                             .set_pdf_end(pdf_end)
                             .set_var_set(var_set));
    wgamma_seq.InsertSequence(wgamma_label, analysis);
  }


  ac::Sequence tp_seq;
  std::string tp_label = "TP";
  if (sequences.count(tp_label)) {
    auto tp_fs = fs.at(tp_label).get();

    if (jsc.count("stitching")) {
      tp_seq.BuildModule(ac::SampleStitching("SampleStitching", jsc["stitching"]));
    }

    tp_seq.BuildModule(ac::EventCounters("EventCounters").set_fs(tp_fs));

    if (is_data) {
      tp_seq.BuildModule(
          ac::LumiMask("LumiMask").set_fs(tp_fs).set_input_file(jsc["data_json"]));
    }

    tp_seq.BuildModule(ac::WGTagAndProbe("WGTagAndProbe")
                             .set_fs(tp_fs)
                             .set_year(jsc["year"])
                             .set_corrections("wgamma/inputs/wgamma_corrections_" + s_year + "_v12.root")
                             .set_is_data(is_data));

    tp_seq.InsertSequence(tp_label, analysis);
  }

  ac::Sequence photp_seq;
  std::string photp_label = "PhotonTP";
  if (sequences.count(photp_label)) {
    auto tp_fs = fs.at(photp_label).get();

    if (jsc.count("stitching")) {
      photp_seq.BuildModule(ac::SampleStitching("SampleStitching", jsc["stitching"]));
    }

    photp_seq.BuildModule(ac::EventCounters("EventCounters").set_fs(tp_fs));

    if (is_data) {
      photp_seq.BuildModule(
          ac::LumiMask("LumiMask").set_fs(tp_fs).set_input_file(jsc["data_json"]));
    }

    photp_seq.BuildModule(ac::WGTagAndProbe("WGTagAndProbe")
                             .set_fs(tp_fs)
                             .set_year(jsc["year"])
                             .set_corrections("wgamma/inputs/wgamma_corrections_" + s_year + "_v11.root")
                             .set_is_data(is_data)
                             .set_do_photons(true));

    photp_seq.InsertSequence(photp_label, analysis);
  }


  ac::Sequence wg_gen_seq;

  if (sequences.count("wg_gen")) {
    if (jsc.count("stitching")) {
      wg_gen_seq.BuildModule(ac::SampleStitching("SampleStitching", jsc["stitching"]));
    }
    wg_gen_seq.BuildModule(ac::WGAnalysis("WGAnalysis")
                              .set_fs(fs.at("wg_gen").get())
                              .set_add_standalone(ac::contains(jsc["attributes"], "add_standalone"))
                              .set_add_rivet(ac::contains(jsc["attributes"], "add_rivet")));
    wg_gen_seq.InsertSequence("wg_gen", analysis);
  }

  ac::Sequence jes_seq;
  if (sequences.count("JES")) {
    auto jes_fs = fs.at("JES").get();
    jes_seq.BuildModule(ac::JESStudy("JESStudy").set_fs(jes_fs));
    jes_seq.InsertSequence("JES", analysis);
  }

  analysis.RunAnalysis();
}
