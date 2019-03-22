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
#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/Photon.h"

namespace ac {

template <typename T, typename Range>
bool contains(Range const& r, T const& value) {
  return std::find(r.begin(), r.end(), value) != r.end();
}

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
bool DescendingTrackPt(Track const* t1, Track const* t2);

// Nicer wrapper for calling a RooFunctor
template <class T>
double RooFunc(T const& func, std::vector<double> const& args) {
  return func->eval(args.data());
}

template <typename T>
std::vector<std::pair<T, T> > MakePairs(std::vector<T> const& collection) {
  unsigned n = collection.size();
  std::vector<std::pair<T, T> > pairVec;
  if (n == 0) return pairVec;
  pairVec.resize((n * (n - 1)) / 2);
  unsigned vecIndex = 0;
  for (unsigned i = 0; i < (n - 1); ++i) {
    for (unsigned j = (i + 1); j < n; ++j, ++vecIndex) {
      pairVec[vecIndex] = (std::pair<T, T>(collection[i], collection[j]));
    }
  }
  return pairVec;
}

template <class T, class U>
std::vector<std::pair<T, U> > MakePairs(std::vector<T> const& collection1,
                                        std::vector<U> const& collection2) {
  unsigned n = collection1.size();
  unsigned m = collection2.size();
  std::vector<std::pair<T, U> > pairVec;
  pairVec.resize(n * m);
  unsigned vecIndex = 0;
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < m; ++j, ++vecIndex) {
      pairVec[vecIndex] = (std::pair<T, U>(collection1[i], collection2[j]));
    }
  }
  return pairVec;
}

// Calculating observables
double DeltaR(ac::Candidate const* c1, ac::Candidate const* c2);
double DeltaRTrack(ac::Track const* c1, ac::Candidate const* c2);
double DeltaRDiTrack(ac::Track const* c1, ac::Track const* c2);
double DeltaRTrackPair(ac::Track const* c1, ac::Track const* c2, ac::Track const* c3);


bool AscendingDR(std::pair<std::pair<unsigned,unsigned>,double> m1, std::pair<std::pair<unsigned,unsigned>,double>m2);
bool DescendingPairPt(std::pair<std::pair<unsigned,unsigned>,double> m1 ,std::pair<std::pair<unsigned,unsigned>,double> m2);

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

// Photon ID/Iso
//  - Option to apply the cut-based photon ID without the charged iso WP or sigmaIetaIeta so we can invert them
//  - year: 2016/2017/2018
//  - wp: 0=loose, 1=medium, 2=tight
bool PhotonIDIso(ac::Photon const* p, unsigned year, unsigned wp, bool apply_charged, bool apply_sigma);

bool ElectronIPCuts(ac::Electron const* e);

// Transverse mass
double MT(Candidate const* cand1, Candidate const* cand2);
}  // namespace ac
#endif
