#ifndef Acorn_Analysis_AnalysisTools_h
#define Acorn_Analysis_AnalysisTools_h
#include <string>
#include <cstdint>
#include <algorithm>
#include "boost/range/algorithm_ext/erase.hpp"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include "TRandom3.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/NTupler/interface/Muon.h"

namespace ac {

// Containers - filtering, sorting etc
template <class Container, class Pred>
Container& keep_if(Container& target, Pred pred) {
  return boost::remove_erase_if(target,
                                [&](typename Container::value_type const& x) { return !pred(x); });
}

template <class Container, class Pred>
Container copy_keep_if(Container const& target, Pred pred) {
  Container res = target;
  keep_if(res, pred);
  return res;
}

bool DescendingPt(Candidate const* c1, Candidate const* c2);

// Nicer wrapper for calling a RooFunctor
template <class T>
double RooFunc(T const& func, std::vector<double> const& args) {
  return func->eval(args.data());
}

// Calculating observables
double DeltaR(ac::Candidate const* c1, ac::Candidate const* c2);

// Trigger matching
bool IsFilterMatchedDR(Candidate const* cand, std::vector<TriggerObject*> const& objs,
                       std::string const& filter, double const& max_dr);

class LookupFilter {
public:
  LookupFilter() = default;
  LookupFilter(std::map<unsigned, std::string> info);

  std::string Lookup(unsigned search) const;

private:
  std::map<unsigned, std::string> info_;
};

// Muon isolation
double MuonPFIso(ac::Muon const* mu);

}  // namespace ac
#endif