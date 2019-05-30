#ifndef Acorn_Analysis_EventCounters_h
#define Acorn_Analysis_EventCounters_h
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

namespace ac {

class EventCounters : public ModuleBase {
 private:
  CLASS_MEMBER(EventCounters, fwlite::TFileService *, fs)

 private:
  TH1D *out;
  std::vector<TH1D *> extra_sets;
  std::vector<bool> extra_sets_relative;

 public:
  EventCounters(std::string const &name);
  EventCounters & AddWeightSet(std::string label, unsigned n_weights, bool is_relative);
  virtual ~EventCounters();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};

}  // namespace ac
#endif
