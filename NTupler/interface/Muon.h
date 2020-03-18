#ifndef Acorn_Muon_h
#define Acorn_Muon_h
#include <map>
#include <string>
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "Acorn/NTupler/interface/Candidate.h"
#include "Rtypes.h"

namespace ac {

/**
 * @brief This class stores a subset of the reco::Muon
 * properties which are most commonly used in analysis.
 *
 * Useful links:
 *  - https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
 */
class Muon : public Candidate {
 public:
  Muon();
  virtual ~Muon();
  virtual void Print() const;

  inline bool isLooseMuon() const { return isLooseMuon_; }
  inline bool isMediumMuon() const { return isMediumMuon_; }
  inline bool isTightMuon() const { return isTightMuon_; }

  inline double dxy() const { return dxy_; }
  inline double dz() const { return dz_; }

  inline float pfIsoSumChargedHadronPt() const { return pfIsoSumChargedHadronPt_; }
  inline float pfIsoSumNeutralHadronEt() const { return pfIsoSumNeutralHadronEt_; }
  inline float pfIsoSumPhotonEt() const { return pfIsoSumPhotonEt_; }
  inline float pfIsoSumPUPt() const { return pfIsoSumPUPt_; }

  inline ROOT::Math::XYZPoint vertex() const { return vertex_; }

  inline void setIsLooseMuon(bool const& isLooseMuon) { isLooseMuon_ = isLooseMuon; }
  inline void setIsMediumMuon(bool const& isMediumMuon) { isMediumMuon_ = isMediumMuon; }
  inline void setIsTightMuon(bool const& isTightMuon) { isTightMuon_ = isTightMuon; }

  inline void setDxy(double const& dxy) { dxy_ = dxy; }
  inline void setDz(double const& dz) { dz_ = dz; }

  inline void setVertex(ROOT::Math::XYZPoint const& vertex) { vertex_ = vertex; }

  inline void SetPfIsoSumChargedHadronPt(float const& pfIsoSumChargedHadronPt) {
    pfIsoSumChargedHadronPt_ = pfIsoSumChargedHadronPt;
  }
  inline void SetPfIsoSumNeutralHadronEt(float const& pfIsoSumNeutralHadronEt) {
    pfIsoSumNeutralHadronEt_ = pfIsoSumNeutralHadronEt;
  }
  inline void SetPfIsoSumPhotonEt(float const& pfIsoSumPhotonEt) {
    pfIsoSumPhotonEt_ = pfIsoSumPhotonEt;
  }
  inline void SetPfIsoSumPUPt(float const& pfIsoSumPUPt) {
    pfIsoSumPUPt_ = pfIsoSumPUPt;
  }

 private:
  bool isLooseMuon_;
  bool isMediumMuon_;
  bool isTightMuon_;

  double dxy_;
  double dz_;

  float pfIsoSumChargedHadronPt_;
  float pfIsoSumNeutralHadronEt_;
  float pfIsoSumPhotonEt_;
  float pfIsoSumPUPt_;

  ROOT::Math::XYZPoint vertex_;
};

}
#endif
