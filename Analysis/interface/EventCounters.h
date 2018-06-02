#ifndef Acorn_Analysis_EventCounters_h
#define Acorn_Analysis_EventCounters_h
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

namespace ac {

class EventCounters : public ModuleBase {
 private:
  CLASS_MEMBER(EventCounters, fwlite::TFileService *, fs)
  TH1D *out;

 public:
  EventCounters(std::string const &name);
  virtual ~EventCounters();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};

}  // namespace ac
#endif
