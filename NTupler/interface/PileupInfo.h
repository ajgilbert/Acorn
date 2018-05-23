#ifndef Acorn_PileupInfo_hh
#define Acorn_PileupInfo_hh
#include <string>
#include <vector>
#include <iostream>
#include "Rtypes.h"

namespace ac {

/**
 * @brief Stores information on the in-time or out-of-time simulated pileup
 * interactions
 */
class PileupInfo {
 public:
  PileupInfo();
  virtual ~PileupInfo();
  virtual void Print() const;


  inline int numInteractions() const { return numInteractions_; }
  inline int bunchCrossing() const { return bunchCrossing_; }
  inline float trueNumInteractions() const { return trueNumInteractions_; }

  inline void setNumInteractions(int const& numInteractions) {
    numInteractions_ = numInteractions;
  }

  inline void setBunchCrossing(int const& bunchCrossing) {
    bunchCrossing_ = bunchCrossing;
  }

  inline void setTrueNumInteractions(float const& trueNumInteractions) {
    trueNumInteractions_ = trueNumInteractions;
  }

 private:
  int numInteractions_;
  int bunchCrossing_;
  float trueNumInteractions_;
};

}
#endif
