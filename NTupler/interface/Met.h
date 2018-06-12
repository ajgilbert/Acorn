#ifndef Acorn_Met_hh
#define Acorn_Met_hh
#include <map>
#include <utility>
#include <vector>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Rtypes.h"

namespace ac {

/**
 * @brief Stores a missing transverse energy object and the corresponding
 * significance and corrections.
 */
class Met : public Candidate {
 public:
  Met();
  virtual ~Met();
  virtual void Print() const;

  inline double sumEt() const { return sumEt_; }

  inline void setSumEt(double const& sumEt) { sumEt_ = sumEt; }

 private:
  double sumEt_;
};

}  // namespace ac
#endif
