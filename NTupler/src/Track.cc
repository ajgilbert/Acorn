#include "../interface/Track.h"

namespace ac {

Track::Track()
    : charge_(0),
      hits_(0),
      pixel_hits_(0),
      quality_(0),
      hits_miss_inner_(0){}

  Track::~Track() {}

  void Track::Print() const {
    std::cout << momentum_ <<std::endl;
    std::cout << ref_point_ <<std::endl;
  }
}
