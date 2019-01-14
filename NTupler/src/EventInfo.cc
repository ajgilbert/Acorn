#include "../interface/EventInfo.h"
#include "boost/format.hpp"

namespace ac {
EventInfo::EventInfo()
    : isRealData_(false),
      event_(0),
      run_(0),
      luminosityBlock_(0),
      bunchCrossing_(0),
      nominalGenWeight_(0.),
      nominalLHEWeight_(0.),
      metfilters_(0),
      numVertices_(0) {}

EventInfo::~EventInfo() {}

void EventInfo::Print(unsigned detail) const {
  using boost::format;
  std::cout << format("%s\n") % std::string(30, '=');
  std::cout << format("%-17s | %10i\n")   % "isRealData"       % isRealData();
  std::cout << format("%-17s | %10i\n")   % "event"            % event();
  std::cout << format("%-17s | %10i\n")   % "luminosityBlock"  % luminosityBlock();
  std::cout << format("%-17s | %10i\n")   % "run"              % run();
  std::cout << format("%-17s | %10i\n")   % "bunchCrossing"    % bunchCrossing();
  std::cout << format("%-17s | %10f\n")   % "nominalGenWeight" % nominalGenWeight();
  std::cout << format("%-17s | %10f\n")   % "nominalLHEWeight" % nominalLHEWeight();
  std::cout << format("%-17s | %10i\n")   % "genWeights[n]"    % genWeights().size();
  if (detail > 0) {
    for (unsigned i = 0; i < genWeights().size(); ++i) {
      std::cout << format("    %-13i | %10f\n") % i % genWeights()[i];
    }
  }
  std::cout << boost::format("%-17s | %10i\n")   % "lheWeights[n]"    % lheWeights().size();
  if (detail > 0) {
    for (auto const& it : lheWeights()) {
      std::cout << format("    %-13i | %10f\n") % it.first % it.second;
    }
  }
//   std::cout << boost::format("%-17s | %10.3f\n") % "jet_rho"        % jet_rho_;
//   std::cout << boost::format("%-17s | %10.3f\n") % "lepton_rho"     % lepton_rho_;
//   std::cout << boost::format("%-17s | %10i\n")   % "good_vertices"  % good_vertices_;
//   std::cout << boost::format("%s\n")      % std::string(30, '-');
//   std::cout << boost::format("%-17s\n")   % "weights";
//   std::cout << boost::format("%s\n")      % std::string(30, '-');
//   SDMap::const_iterator it = weights_.begin();
//   SBMap::const_iterator its = weight_status_.begin();
//   for (; it != weights_.end() && its != weight_status_.end(); ++it, ++its) {
//     std::cout << boost::format("%-17s | %6.3f %3i\n") % it->first % it->second %
//                      its->second;
//   }
//   if (filters_.size() > 0) {
//     std::cout << boost::format("%s\n") % std::string(30, '-');
//     std::cout << boost::format("%-17s\n")   % "filters";
//     std::cout << boost::format("%s\n")      % std::string(30, '-');
//     TBMap::const_iterator itf = filters_.begin();
//     for (; itf != filters_.end(); ++itf) {
//       std::cout << boost::format("%-21s | %6.3f\n") % itf->first % itf->second;
//     }
//   }
//   std::cout << boost::format("%s\n")      % std::string(30, '=');
}
}
