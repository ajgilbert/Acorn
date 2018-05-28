#include "../interface/GenParticle.h"
#include "../interface/Reduction.h"
#include "boost/format.hpp"
namespace ac {
// Constructors/Destructors
GenParticle::GenParticle() : index_(0), pdgId_(0), status_(0), spin_(0.) {}

GenParticle::~GenParticle() {}

void GenParticle::Print() const {
  std::cout << (boost::format("idx: %-4i  st: %-3i  id: %4i  %-40s  M: %f  sp: %f  flags: %s\n") %
                this->index() % this->status() % this->pdgId()) %
                   this->vector() % this->M() % this->spin() % toBinaryString(this->statusFlags_);
  if (this->mothers().size()) {
    std::cout << "  mothers:  ";
    for (unsigned i = 0; i < this->mothers().size(); ++i) {
      std::cout << " " << this->mothers().at(i);
    }
    std::cout << "\n";
  }
  if (this->daughters().size()) {
    std::cout << "  daughters:";
    for (unsigned i = 0; i < this->daughters().size(); ++i) {
      std::cout << " " << this->daughters().at(i);
    }
    std::cout << "\n";
  }
}
}
