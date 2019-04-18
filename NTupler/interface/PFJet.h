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

  inline void setPassesJetID(bool const& passesJetID) { passesJetID_ = passesJetID; }

 private:
  bool passesJetID_;

};

}
#endif
