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

  inline double scEta() const { return scEta_; }
  inline float hadTowOverEm() const { return hadTowOverEm_; }
  inline float full5x5SigmaIetaIeta() const { return full5x5SigmaIetaIeta_; }

  inline double chargedIso() const { return chargedIso_; }
  inline double neutralHadronIso() const { return neutralHadronIso_; }
  inline double photonIso() const { return photonIso_; }

  inline std::vector<float> const& energyCorrections() const { return energyCorrections_; }


  inline void setIsLooseIdPhoton(bool const& isLooseIdPhoton) { isLooseIdPhoton_ = isLooseIdPhoton; }
  inline void setIsMediumIdPhoton(bool const& isMediumIdPhoton) { isMediumIdPhoton_ = isMediumIdPhoton; }
  inline void setIsTightIdPhoton(bool const& isTightIdPhoton) { isTightIdPhoton_ = isTightIdPhoton; }

  inline void setPassElectronVeto(bool const& passElectronVeto) { passElectronVeto_ = passElectronVeto; }
  inline void setHasPixelSeed(bool const& hasPixelSeed) { hasPixelSeed_ = hasPixelSeed; }

  inline void setScEta(double const& scEta) { scEta_ = scEta; }
  inline void setHadTowOverEm(float const& hadTowOverEm) { hadTowOverEm_ = hadTowOverEm; }
  inline void setFull5x5SigmaIetaIeta(float const& full5x5SigmaIetaIeta) { full5x5SigmaIetaIeta_ = full5x5SigmaIetaIeta; }

  inline void setChargedIso(double const& chargedIso) { chargedIso_ = chargedIso; }
  inline void setNeutralHadronIso(double const& neutralHadronIso) { neutralHadronIso_ = neutralHadronIso; }
  inline void setPhotonIso(double const& photonIso) { photonIso_ = photonIso; }

  inline void setEnergyCorrections(std::vector<float> const& energyCorrections) { energyCorrections_ = energyCorrections; }


 private:
  bool isLooseIdPhoton_;
  bool isMediumIdPhoton_;
  bool isTightIdPhoton_;

  bool passElectronVeto_;
  bool hasPixelSeed_;

  double scEta_;
  float hadTowOverEm_;
  float full5x5SigmaIetaIeta_;

  double chargedIso_;
  double neutralHadronIso_;
  double photonIso_;

  std::vector<float> energyCorrections_;
};
}
#endif
