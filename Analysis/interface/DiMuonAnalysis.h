#ifndef Acorn_Analysis_DiMuonAnalysis_h
#define Acorn_Analysis_DiMuonAnalysis_h
#include <string>
#include <cstdint>
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

namespace ac {

template <class Container, class Pred>
Container& keep_if(Container& target, Pred pred) {
  return boost::remove_erase_if(target,
                                [&](typename Container::value_type const& x) { return !pred(x); });
}

/**
 * Copy a container, filter elements, then return this copy
 */
template <class Container, class Pred>
Container copy_keep_if(Container const& target, Pred pred) {
  Container res = target;
  keep_if(res, pred);
  return res;
}

class DiMuonAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(DiMuonAnalysis, fwlite::TFileService*, fs)

  TTree* tree_;
  float pt_1_;
  float pt_2_;
  float m_ll_;

 public:
  DiMuonAnalysis(std::string const& name);
  virtual ~DiMuonAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
