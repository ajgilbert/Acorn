#include "../interface/RunLumiMap.h"
#include <iostream>
#include "boost/lexical_cast.hpp"
#include "Acorn/NTupler/interface/json.hpp"

RunLumiMap::RunLumiMap() {}

RunLumiMap::RunLumiMap(const char *name, const char *title) : TNamed(name, title) {}

RunLumiMap::~RunLumiMap() {}

Long64_t RunLumiMap::Merge(TCollection* coll) {
  for (auto it = coll->begin(); it != coll->end(); ++it) {
    RunLumiMap const* other = dynamic_cast<RunLumiMap const*>(*it);
    if (other) {
      AppendMap(other->lsmap_);
    }
  }
  return 0;
}

void RunLumiMap::Remove(unsigned run, unsigned ls) {
  auto it = lsmap_.find(run);
  if (it == lsmap_.end()) {
    return;
  } else {
    auto it2 = it->second.find(ls);
    if (it2 == it->second.end()) {
      return;
    } else {
      it->second.erase(it2);
      if (it->second.size() == 0) {
        lsmap_.erase(it);
      }
    }
  }
}

void RunLumiMap::Add(unsigned run, unsigned ls) {
  lsmap_[run].insert(ls);
}

bool RunLumiMap::InMap(unsigned run, unsigned ls) const {
  auto it = lsmap_.find(run);
  if (it == lsmap_.end()) {
    return false;
  } else {
    auto it2 = it->second.find(ls);
    if (it2 == it->second.end()) {
      return false;
    }
  }
  return true;
}

std::map<unsigned, std::set<unsigned>> const& RunLumiMap::GetMap() const {
  return lsmap_;
}
void RunLumiMap::ClearMap() {
  lsmap_.clear();
}

void RunLumiMap::ResetMap(std::map<unsigned, std::set<unsigned>> const& newmap) {
  ClearMap();
  lsmap_ = newmap;
}
void RunLumiMap::AppendMap(std::map<unsigned, std::set<unsigned>> const& newmap) {
  for (auto rit = newmap.begin(); rit != newmap.end(); ++rit) {
    for (auto lit = rit->second.begin(); lit != rit->second.end(); ++lit) {
      this->Add(rit->first, *lit);
    }
  }
}

void RunLumiMap::PrintMap() const {
  for (auto rit = lsmap_.begin(); rit != lsmap_.end(); ++rit) {
    std::cout << rit->first << ": [";
    for (auto lit = rit->second.begin(); lit != rit->second.end(); ++lit) {
      std::cout << *lit;
      auto litcopy = lit;
      if (++litcopy != rit->second.end()) {
        std::cout << ", ";
      }
    }
    std::cout << "]\n";
  }
}

void RunLumiMap::AppendFromJsonString(std::string const& js_str) {
  nlohmann::json js = nlohmann::json::parse(js_str);
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
          Add(run, i);
      }
    }
  }
}

std::string RunLumiMap::AsJsonString() const {
  nlohmann::json js;
  for (auto const& info : GetMap()) {
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
  return js.dump();
}

ClassImp(RunLumiMap)
