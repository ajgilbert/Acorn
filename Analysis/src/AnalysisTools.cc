#include "Acorn/Analysis/interface/AnalysisTools.h"
#include <algorithm>
#include <map>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/Track.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "Acorn/NTupler/interface/TriggerObject.h"
#include "Acorn/NTupler/interface/city.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMath.h"
#include "boost/lexical_cast.hpp"

namespace ac {

// Containers - filtering, sorting etc
bool DescendingPt(Candidate const* c1, Candidate const* c2) { return c1->pt() > c2->pt(); }
bool DescendingTrackPt(Track const* t1, Track const* t2) { return t1->pt() > t2->pt(); }

bool AscendingDR(std::pair<std::pair<unsigned,unsigned>,double> m1, std::pair<std::pair<unsigned,unsigned>,double>m2) { return m1.second < m2.second; }
bool DescendingPairPt(std::pair<std::pair<unsigned,unsigned>,double> m1, std::pair<std::pair<unsigned,unsigned>,double> m2) {return m1.second > m2.second;}


// Calculating observables
double DeltaR(ac::Candidate const* c1, ac::Candidate const* c2) {
  return ROOT::Math::VectorUtil::DeltaR(c1->vector(), c2->vector());
}

double DeltaRTrack(ac::Track const* c1, ac::Candidate const* c2) {
  return ROOT::Math::VectorUtil::DeltaR(c1->vector(), c2->vector());
}

double DeltaRDiTrack(ac::Track const* c1, ac::Track const* c2) {
  return ROOT::Math::VectorUtil::DeltaR(c1->vector(), c2->vector());
}

double DeltaRTrackPair(ac::Track const* c1, ac::Track const* c2, ac::Track const* c3){
 return ROOT::Math::VectorUtil::DeltaR(c1->vector()+c2->vector(),c3->vector());
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
    std::cerr << "Error in LookupFilter::Lookup, map is empty\n";
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

double MTcluster(Candidate const* cand1a, Candidate const* cand1b, Candidate const* cand2) {
  auto cand1 = cand1a->vector() + cand1b->vector();
  auto full_system = cand1 + cand2->vector();
  double mTcluster2 = std::pow(std::sqrt(cand1.M2() + cand1.perp2()) + cand2->pt(), 2) - full_system.perp2();
  if (mTcluster2 > 0) {
    return std::sqrt(mTcluster2);
  } else {
    std::cerr << "Cluster transverse mass would be negative! Returning 0.0" << std::endl;
  }
  return 0.0;
}


bool ElectronIsoFall17V2(ac::Electron const* e, unsigned wp) {
  bool eb = std::abs(e->scEta()) < 1.479;
  double iso = e->relativeEAIso();
  double pt = e->pt();
  if (eb) {
    if (wp == 0) {  // Veto
      return iso < 0.198 + 0.506 / pt;
    }
    if (wp == 1) {  // Loose
      return iso < 0.112 + 0.506 / pt;
    }
    if (wp == 2) {  // Medium
      return iso < 0.0478 + 0.506 / pt;
    }
    if (wp == 3) {  // Tight
      return iso < 0.0287 + 0.506 / pt;
    }
  } else {
    if (wp == 0) {  // Veto
      return iso < 0.203 + 0.963 / pt;
    }
    if (wp == 1) {  // Loose
      return iso < 0.108 + 0.963 / pt;
    }
    if (wp == 2) {  // Medium
      return iso < 0.0658 + 0.963 / pt;
    }
    if (wp == 3) {  // Tight
      return iso < 0.0445 + 0.963 / pt;
    }
  }
  return true;
}

bool PhotonIDIso(ac::Photon const* p, unsigned year, unsigned wp, bool apply_charged, bool apply_sigma) {
  double pt = p->pt();
  bool eb = std::abs(p->scEta()) < 1.4442;
  double hovere = 0.;
  double sigma = 0.;
  double charged = 0.;
  double neutral = 0.;
  double photon = 0.;
  if (year == 2016 || year == 2017 || year == 2018) {
    if (eb && wp == 0) {
      hovere = 0.04596;
      sigma = 0.0106;
      charged = 1.694;
      neutral = 24.032 + 0.01512 * pt + 2.259e-05 * pt * pt;
      photon = 2.876 + 0.004017 * pt;
    }
    if (eb && wp == 1) {
      hovere = 0.02197;
      sigma = 0.01015;
      charged = 1.141;
      neutral = 1.189 + 0.01512 * pt + 2.259e-05 * pt * pt;
      photon = 2.08 + 0.004017 * pt;
    }
    if (eb && wp == 2) {
      hovere = 0.02148;
      sigma = 0.00996;
      charged = 0.65;
      neutral = 0.317 + 0.01512 * pt + 2.259e-05 * pt * pt;
      photon = 2.044 + 0.004017 * pt;
    }
    if (!eb && wp == 0) {
      hovere = 0.0590;
      sigma = 0.0272;
      charged = 2.089;
      neutral = 19.722 + 0.0117 * pt + 2.3e-05 * pt * pt;
      photon = 4.162 +  0.0037 * pt;
    }
    if (!eb && wp == 1) {
      hovere = 0.0326;
      sigma = 0.0272;
      charged = 1.051;
      neutral = 1.715 + 0.0117 * pt + 2.3e-05 * pt * pt;
      photon = 3.863 +  0.0037 * pt;
    }
    if (!eb && wp == 2) {
      hovere = 0.0321;
      sigma = 0.0271;
      charged = 0.517;
      neutral = 0.586 + 0.0117 * pt + 2.3e-05 * pt * pt;
      photon = 2.617 +  0.0037 * pt;
    }
  }
  return p->hadTowOverEm() < hovere && (apply_sigma ? (p->full5x5SigmaIetaIeta() < sigma) : true) &&
         (apply_charged ? (p->chargedIso() < charged) : true) && p->neutralHadronIso() < neutral &&
         p->photonIso() < photon;
}

bool ElectronIPCuts(ac::Electron const* e) {
  return (fabs(e->eta()) < 1.4442 && fabs(e->dxy()) < 0.05 && fabs(e->dz()) < 0.1) ||
         (fabs(e->eta()) >= 1.4442 && fabs(e->dxy()) < 0.10 && fabs(e->dz()) < 0.2);
}


}  // namespace ac
