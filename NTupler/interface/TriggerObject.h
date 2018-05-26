#ifndef Acorn_TriggerObject_h
#define Acorn_TriggerObject_h
#include <vector>
#include "Acorn/NTupler/interface/Candidate.h"
#include "Rtypes.h"

namespace ac {

/**
 * @brief Stores the four-momentum of a trigger object as well as a list of the
 * (hashed) filter labels the object was accepted by
 */
class TriggerObject : public Candidate {
 private:
 public:
  TriggerObject();
  virtual ~TriggerObject();
  virtual void Print() const;

  inline std::vector<std::uint64_t> const& filters() const { return filters_; }

  inline void setFilters(std::vector<std::uint64_t> const& filters) {
    filters_ = filters;
  }

 private:
  std::vector<std::uint64_t> filters_;
};

}
#endif
