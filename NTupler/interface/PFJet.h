#ifndef Acorn_PFJet_h
#define Acorn_PFJet_h
#include <string>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Rtypes.h"

namespace ac {

class PFJet : public Candidate {
 public:
  PFJet();
  virtual ~PFJet();
  virtual void Print() const;

  inline bool passesJetID() const { return passesJetID_; }
  inline float deepCSVDiscriminatorBvsAll() const { return deepCSVDiscriminatorBvsAll_; }

  inline void setPassesJetID(bool const& passesJetID) { passesJetID_ = passesJetID; }
  inline void setDeepCSVDiscriminatorBvsAll(float const& deepCSVDiscriminatorBvsAll) { deepCSVDiscriminatorBvsAll_ = deepCSVDiscriminatorBvsAll; }

 private:
  bool passesJetID_;
  float deepCSVDiscriminatorBvsAll_;

};

}
#endif
