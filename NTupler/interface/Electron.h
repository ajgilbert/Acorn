#ifndef Acorn_Electron_h
#define Acorn_Electron_h
#include <map>
#include <string>
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "Acorn/NTupler/interface/Candidate.h"
#include "Rtypes.h"

namespace ac {
/**
 * @brief This class stores a subset of the reco::Electron
 * properties which are most commonly used in analysis.
 *
 * Useful links:
 */


class Electron : public Candidate {
 public:
  Electron();
  virtual ~Electron();
  virtual void Print() const;

  inline bool isCutBasedVetoElectron() const { return isCutBasedVetoElectron_; }
  inline bool isCutBasedLooseElectron() const { return isCutBasedLooseElectron_; }
  inline bool isCutBasedMediumElectron() const { return isCutBasedMediumElectron_; }
  inline bool isCutBasedTightElectron() const { return isCutBasedTightElectron_; }

  inline bool isMVAwp90Electron() const { return isMVAwp90Electron_; }
  inline bool isMVAwp80Electron() const { return isMVAwp80Electron_; }

  inline bool isHEEPElectron() const { return isHEEPElectron_; }

  inline double relativeEAIso() const { return relativeEAIso_; }

  inline double dxy() const { return dxy_; }
  inline double dz() const { return dz_; }

  inline ROOT::Math::XYZPoint vertex() const { return vertex_; }

  inline void setIsCutBasedVetoElectron(bool const& isCutBasedVetoElectron) { isCutBasedVetoElectron_ = isCutBasedVetoElectron; }
  inline void setIsCutBasedLooseElectron(bool const& isCutBasedLooseElectron) { isCutBasedLooseElectron_ = isCutBasedLooseElectron; }
  inline void setIsCutBasedMediumElectron(bool const& isCutBasedMediumElectron) { isCutBasedMediumElectron_ = isCutBasedMediumElectron; }
  inline void setIsCutBasedTightElectron(bool const& isCutBasedTightElectron) { isCutBasedTightElectron_ = isCutBasedTightElectron; }

  inline void setIsMVAwp90Electron(bool const& isMVAwp90Electron) { isMVAwp90Electron_ = isMVAwp90Electron; }
  inline void setIsMVAwp80Electron(bool const& isMVAwp80Electron) { isMVAwp80Electron_ = isMVAwp80Electron; }

  inline void setIsHEEPElectron(bool const& isHEEPElectron) { isHEEPElectron_ = isHEEPElectron; }

  inline void setRelativeEAIso(double const& relativeEAIso) { relativeEAIso_ = relativeEAIso; }

  inline void setDxy(double const& dxy) { dxy_ = dxy; }
  inline void setDz(double const& dz) { dz_ = dz; }

  inline void setVertex(ROOT::Math::XYZPoint const& vertex) { vertex_ = vertex; }


 private:
  bool isCutBasedVetoElectron_;
  bool isCutBasedLooseElectron_;
  bool isCutBasedMediumElectron_;
  bool isCutBasedTightElectron_;

  bool isMVAwp90Electron_;
  bool isMVAwp80Electron_;

  bool isHEEPElectron_;

  double relativeEAIso_;

  double dxy_;
  double dz_;

  ROOT::Math::XYZPoint vertex_;

};

}
#endif
