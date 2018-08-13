#include "Acorn/Analysis/interface/AnalysisTools.h"
#include <algorithm>
#include <map>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/NTupler/interface/city.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMath.h"
#include "boost/lexical_cast.hpp"

namespace ac {

// Containers - filtering, sorting etc
bool DescendingPt(Candidate const* c1, Candidate const* c2) { return c1->pt() > c2->pt(); }

// Calculating observables
double DeltaR(ac::Candidate const* c1, ac::Candidate const* c2) {
  return ROOT::Math::VectorUtil::DeltaR(c1->vector(), c2->vector());
}

bool IsFilterMatchedDR(Candidate const* cand, std::vector<TriggerObject*> const& objs,
                     std::string const& filter, double const& max_dr) {
  std::size_t hash = CityHash64(filter);
  for (unsigned i = 0; i < objs.size(); ++i) {
    std::vector<std::uint64_t> const& labels = objs[i]->filters();
    if (std::find(labels.begin(), labels.end(), hash) == labels.end()) continue;
    if (DeltaR(cand, objs[i]) < max_dr) return true;
  }
  return false;
}

LookupFilter::LookupFilter(std::map<unsigned, std::string> info) {
  info_ = info;
}
std::string LookupFilter::Lookup(unsigned search) const {
  if (info_.size() == 0) {
    std::cerr << ">> Error in LookupFilter::Lookup, map is empty\n";
    return "";
  } else {
    // An iterator to the the first element in the container whose key is considered to go after k,
    // or map::end if no keys are considered to go after k.
    // I.e. if searching for 0 in { 1, 8, 20 } -> returns it to 1 -> OK
    //                       1 in { 1, 8, 20 } -> returns it to 2 -> decreate by 1
    //                       4 in { 1, 8, 20 } -> returns it end() -> decrease by 1

    auto it = info_.upper_bound(search);
    if (it != info_.begin()) {
      --it;
    }
    return it->second;
  }
}

// Muon isolation
double MuonPFIso(ac::Muon const* m) {
  return (m->pfIsoSumChargedHadronPt() +
          std::max(
              0., m->pfIsoSumNeutralHadronEt() + m->pfIsoSumPhotonEt() - 0.5 * m->pfIsoSumPUPt())) /
         m->pt();
}

double MT(Candidate const* cand1, Candidate const* cand2) {
  double mt = 2. * cand1->pt() * cand2->pt() *
              (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(cand1->vector(), cand2->vector())));
  if (mt > 0) {
    return std::sqrt(mt);
  } else {
    std::cerr << "Transverse mass would be negative! Returning 0.0" << std::endl;
  }
  return 0.0;
}

bool PhotonIDIso(ac::Photon const* p, unsigned year, unsigned wp, bool apply_charged) {
  double pt = p->pt();
  bool eb = std::abs(p->scEta()) < 1.4442;
  double hovere = 0.;
  double sigma = 0.;
  double charged = 0.;
  double neutral = 0.;
  double photon = 0.;
  if (year == 2016) {
    if (eb && wp == 0) {
      hovere = 0.0597;
      sigma = 0.01031;
      charged = 1.295;
      neutral = 10.910 + 0.0148 * pt + 0.000017 * pt * pt;
      photon = 3.630 + 0.0047 * pt;
    }
    if (eb && wp == 1) {
      hovere = 0.0396;
      sigma = 0.01022;
      charged = 0.441;
      neutral = 2.725 + 0.0148 * pt + 0.000017 * pt * pt;
      photon = 2.571 + 0.0047 * pt;
    }
    if (eb && wp == 2) {
      hovere = 0.0269;
      sigma = 0.00994;
      charged = 0.202;
      neutral = 0.264 + 0.0148 * pt + 0.000017 * pt * pt;
      photon = 2.362 + 0.0047 * pt;
    }
    if (!eb && wp == 0) {
      hovere = 0.0481;
      sigma = 0.03013;
      charged = 1.011;
      neutral = 5.931 + 0.0163 * pt + 0.000014 * pt * pt;
      photon = 6.641 + 0.0034 * pt;
    }
    if (!eb && wp == 1) {
      hovere = 0.0219;
      sigma = 0.03001;
      charged = 0.442;
      neutral = 1.715 + 0.0163 * pt + 0.000014 * pt * pt;
      photon = 3.863 + 0.0034 * pt;
    }
    if (!eb && wp == 2) {
      hovere = 0.0213;
      sigma = 0.03000;
      charged = 0.034;
      neutral = 0.586 + 0.0163 * pt + 0.000014 * pt * pt;
      photon = 2.617 + 0.0034 * pt;
    }
  }
  return p->hadTowOverEm() < hovere && p->full5x5SigmaIetaIeta() < sigma &&
         (apply_charged ? (p->chargedIso() < charged) : true) && p->neutralHadronIso() < neutral &&
         p->photonIso() < photon;
}

}  // namespace ac
