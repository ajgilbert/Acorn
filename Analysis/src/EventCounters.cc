#include "Acorn/Analysis/interface/EventCounters.h"
#include "Acorn/NTupler/interface/EventInfo.h"

namespace ac {
EventCounters::EventCounters(std::string const& name) : ModuleBase(name) { fs_ = NULL; }

EventCounters::~EventCounters() { ; }

int EventCounters::PreAnalysis() {
  std::vector<std::string> labels = {"unweighted", "weighted", "pos", "neg"};
  unsigned n = labels.size();
  out = fs_->make<TH1D>("counters", "counters", n, 0., double(n));
  for (unsigned i = 0; i < n; ++i) {
    out->GetXaxis()->SetBinLabel(i+1, labels[i].c_str());
  }
  return 0;
}

int EventCounters::Execute(TreeEvent* event) {
  EventInfo const* info = event->GetPtr<EventInfo>("eventInfo");
  // unweighted
  out->Fill(0.5, 1.0);
  // weighted
  // double w = info->totalWeight();
  double w = info->nominalGenWeight() >= 0. ? +1. : -1.;
  out->Fill(1.5, w);
  if (w >= 0.) {
    // pos
    out->Fill(2.5, w);
  } else {
    // neg
    out->Fill(3.5, w);
  }
  return 0;
}

int EventCounters::PostAnalysis() { return 0; }

void EventCounters::PrintInfo() { ; }
}  // namespace ac
