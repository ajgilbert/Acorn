#ifndef Acorn_Photon_h
#define Acorn_Photon_h
#include <vector>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Rtypes.h"

namespace ac {

/**
 * @brief This class stores a subset of the reco::Photon
 * properties which are most commonly used in analysis.
 */
class Photon : public Candidate {
 public:
  Photon();
  virtual ~Photon();
  virtual void Print() const;

  inline bool isLooseIdPhoton() const { return isLooseIdPhoton_; }
  inline bool isMediumIdPhoton() const { return isMediumIdPhoton_; }
  inline bool isTightIdPhoton() const { return isTightIdPhoton_; }

  inline bool passElectronVeto() const { return passElectronVeto_; }
  inline bool hasPixelSeed() const { return hasPixelSeed_; }

  inline void setIsLooseIdPhoton(bool const& isLooseIdPhoton) { isLooseIdPhoton_ = isLooseIdPhoton; }
  inline void setIsMediumIdPhoton(bool const& isMediumIdPhoton) { isMediumIdPhoton_ = isMediumIdPhoton; }
  inline void setIsTightIdPhoton(bool const& isTightIdPhoton) { isTightIdPhoton_ = isTightIdPhoton; }

  inline void setPassElectronVeto(bool const& passElectronVeto) { passElectronVeto_ = passElectronVeto; }
  inline void setHasPixelSeed(bool const& hasPixelSeed) { hasPixelSeed_ = hasPixelSeed; }

 private:
  bool isLooseIdPhoton_;
  bool isMediumIdPhoton_;
  bool isTightIdPhoton_;

  bool passElectronVeto_;
  bool hasPixelSeed_;
};
}
#endif