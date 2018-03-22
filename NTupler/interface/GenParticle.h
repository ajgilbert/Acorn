#ifndef Acorn_GenParticle_hh
#define Acorn_GenParticle_hh
#include <vector>
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "Acorn/NTupler/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "Rtypes.h"


namespace ac {

/**
 * @brief Stores the basic properties of generator-level particles as well as
 * mother-daughter relations with other particles
 */
class GenParticle : public Candidate {
 public:
  GenParticle();
  virtual ~GenParticle();
  virtual void Print() const;

  /// @name Properties
  /**@{*/
  /// The index position of the particle in the original list
  inline int index() const { return index_; }

  /// PDG number to identify the particle type, see
  /// [this link] (http://pdg.lbl.gov/2002/montecarlorpp.pdf)
  inline int pdgid() const { return pdgid_; }

  /// The generator-dependent particle status
  inline int status() const { return status_; }

  /// A vector of ic::GenParticle::index() values that identify the mother
  /// particles
  inline std::vector<int> const& mothers() const { return mothers_; }

  /// A vector of ic::GenParticle::index() values that identify the daughter
  /// particles
  inline std::vector<int> const& daughters() const { return daughters_; }

  /// A genstatusflags object to give information about the production of the particle
  inline reco::GenStatusFlags statusFlags() const { reco::GenStatusFlags flags; flags.flags_ = std::bitset<15>(statusFlags_); return flags; }
  /**@}*/

  /// @name Setters
  /**@{*/
  /// @copybrief index()
  inline void set_index(int const& index) { index_ = index; }

  /// @copybrief pdgid()
  inline void set_pdgid(int const& pdgid) { pdgid_ = pdgid; }

  /// @copybrief status()
  inline void set_status(int const& status) { status_ = status; }

  /// @copybrief mothers()
  inline void set_mothers(std::vector<int> const& mothers) {
    mothers_ = mothers;
  }

  /// @copybrief daughters()
  inline void set_daughters(std::vector<int> const& daughters) {
    daughters_ = daughters;
  }

  /// @copybrief statusFlags()
  inline void set_statusFlags(reco::GenStatusFlags const& statusFlags) {
    statusFlags_ = statusFlags.flags_.to_ulong();
  }
  /**@}*/

 private:
  int index_;
  int pdgid_;
  int status_;
  std::vector<int> mothers_;
  std::vector<int> daughters_;
  int16_t statusFlags_; // results in a 30% space saving compared to bitset<15>!!
};
}
#endif
