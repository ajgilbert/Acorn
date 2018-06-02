#include "Acorn/Analysis/interface/LumiMask.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "boost/lexical_cast.hpp"
#include <fstream>

namespace ac {

LumiMask::LumiMask(std::string const& name) : ModuleBase(name) {
  fs_ = nullptr;
}

LumiMask::~LumiMask() { ; }

int LumiMask::PreAnalysis() {
  PrintHeader("LumiMask");
  if (input_file_ != "") {
    PrintArg("input_file", input_file_);
    nlohmann::json js;
    std::fstream input;
    input.open(input_file_);
    input >> js;
    FillRunLumiMapFromJson(input_json_, js);
  } else {
    PrintArg("input_file", "-");
  }
  accept_json_ = fs_->make<RunLumiMap>("accept_json", "accept_json");
  reject_json_ = fs_->make<RunLumiMap>("reject_json", "reject_json");
  all_json_ = fs_->make<RunLumiMap>("all_json", "all_json");

  return 0;
}

int LumiMask::Execute(TreeEvent* event) {
  EventInfo const* info = event->GetPtr<EventInfo>("eventInfo");
  unsigned run = info->run();
  unsigned ls = info->luminosityBlock();
  all_json_->Add(run, ls);
  bool accept = true;
  if (input_file_ != "") {
    if (!input_json_.InMap(run, ls)) {
      accept = false;
    }
  }
  if (accept) {
    accept_json_->Add(run, ls);
    return 0;
  } else {
    reject_json_->Add(run, ls);
    return 1;
  }

  return 0;
}
int LumiMask::PostAnalysis() {
  return 0;
}

void FillRunLumiMapFromJson(RunLumiMap & rlmap, nlohmann::json const& js) {
  for (auto const& item : js.items()) {
    nlohmann::json const& run_js = item.value();
    unsigned run = boost::lexical_cast<unsigned>(item.key());
    for (unsigned i = 0; i < run_js.size(); ++i) {
      if (run_js[i].size() != 2) {
        throw std::runtime_error(
            "[LumiMask] Lumi range not in the form [X,Y]");
      }
      unsigned range_min = run_js[i][0];
      unsigned range_max = run_js[i][1];
      if (range_max < range_min) {
        throw std::runtime_error(
            "[LumiMask] Have lumi range [X,Y] where Y < X");
      }
      for (unsigned i = range_min; i <= range_max; ++i) {
          rlmap.Add(run, i);
      }
    }
  }
}

nlohmann::json JsonFromRunLumiMap(RunLumiMap const& rlmap) {
  nlohmann::json js;
  for (auto const& info : rlmap.GetMap()) {
    auto const& run = boost::lexical_cast<std::string>(info.first);
    auto const& lumis = info.second;
    if (lumis.size() == 0) continue;
    auto it = lumis.begin();
    unsigned first_val = *it;
    unsigned last_val = *it;
    while (it != lumis.end()) {
      ++it;
      if (it != lumis.end() && *it == last_val + 1) {
        last_val = *it;
      } else {
        nlohmann::json lumi_pair;
        lumi_pair.push_back(first_val);
        lumi_pair.push_back(last_val);
        js[run].push_back(lumi_pair);
        if (it != lumis.end()) {
          first_val = *it;
          last_val = *it;
        }
      }
    }
  }
  return js;
}

void LumiMask::PrintInfo() { ; }
}
