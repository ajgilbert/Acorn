#ifndef Acorn_Candidate_hh
#define Acorn_Candidate_hh
#include <vector>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Rtypes.h"

namespace ac {

/**
 * @brief Stores a four-momentum, charge and identifier, and is the base class
 * for most other physics objects. */
class Candidate {
 private:
  typedef ROOT::Math::PtEtaPhiMVector Vector;

 public:
  Candidate();
  virtual ~Candidate();
  virtual void Print() const;

  /// @name Properties
  /**@{*/
  /// Four-momentum
  /** Returns the underlying ROOT
   * [PtEtaPhiEVector](http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/classROOT_1_1Math_1_1LorentzVector.html),
   * from which a large number of kinematic variables can be calculated.
   * Please note, specific getter functions are provided by the Candidate
   * class for the most commonly used variables: pt, eta, phi, energy and
   * mass. */
  inline Vector const& vector() const { return vector_; }

  /// Unique identifier
  /** Each Candidate object (and any object deriving from it) typically has
   * a unique id assigned by its Producer. This provides a lightweight
   * mechanism for objects to refer to each other. For example, A GenJet may
   * store the id of each GenParticle constituents. As another example, a
   * PFJet may store the id of its constituent tracks. */
  inline std::size_t id() const { return id_; }

  /// Direct access to the \f$ p_{\text{T}} \f$
  inline double pt() const { return vector_.Pt(); }

  /// Direct access to \f$ \eta \f$
  inline double eta() const { return vector_.Eta(); }

  /// Direct access to \f$ \phi \f$
  inline double phi() const { return vector_.Phi(); }

  /// Direct access to the energy
  inline double energy() const { return vector_.E(); }

  /// Electric charge
  inline int charge() const { return charge_; }

  /// Calculates the candidate mass
  inline double M() const { return vector_.M(); }
  /**@}*/

  /// @name Setters
  /**@{*/
  /// @copybrief vector()
  inline void setVector(Vector const& vector) { vector_ = vector; }

  /// @copybrief id()
  inline void setId(std::size_t const& id) { id_ = id; }

  /// @copybrief pt()
  inline void setPt(double const& pt) { vector_.SetPt(pt); }

  /// @copybrief eta()
  inline void setEta(double const& eta) { vector_.SetEta(eta); }

  /// @copybrief phi()
  inline void setPhi(double const& phi) { vector_.SetPhi(phi); }

  /// @copybrief energy()
  inline void setM(double const& mass) { vector_.SetM(mass); }

  /// @copybrief charge()
  inline void setCharge(int const& charge) { charge_ = charge; }
  /**@}*/

 private:
  Vector vector_;
  std::size_t id_;
  int charge_;
};

}
/** \example plugins/ICCandidateProducer.hh */
#endif
