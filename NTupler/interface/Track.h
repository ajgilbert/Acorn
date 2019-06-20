#ifndef Acorn_Track_h
#define Acorn_Track_h
#include <map>
#include <string>
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Rtypes.h"

namespace ac {
/**
 * @brief This class stores a subset of the reco::Track
 * properties.
 *
 */

class Track {
 private:
  typedef ROOT::Math::RhoEtaPhiVector ThreeVector;
  typedef ROOT::Math::PtEtaPhiEVector Vector;
  typedef ROOT::Math::XYZPoint Point;


 public:
  Track();
  virtual ~Track();
  virtual void Print() const;

  inline ThreeVector const& momentum() const { return momentum_; }
  inline Vector vector(double mass=0.13957018) const {
    return Vector(ROOT::Math::PtEtaPhiMVector(pt(), eta(), phi(), mass));
  }

  inline Point const& ref_point() const { return ref_point_; }

  inline double pt() const { return momentum_.Rho(); }
  inline double energy() const {return momentum_.r(); }
  inline double eta() const { return momentum_.Eta(); }
  inline double phi() const { return momentum_.Phi(); }
  inline double vx() const { return ref_point_.x(); } 
  inline double vy() const { return ref_point_.y(); }
  inline double vz() const { return ref_point_.z(); }
  inline int hits() const { return pixel_hits_;}
  inline double dxy(Point const& point) const {
    return (-(vx() - point.x()) * momentum().y() +
            (vy() - point.y()) * momentum().z()) /
           pt();
  }
  inline double dz(Point const& point) const {
    return (vz() - point.z()) -
           ((vx() - point.x()) * momentum().x() +
            (vy() - point.y()) * momentum().y()) /
               pt() * momentum().z() / pt();
  }
  inline int charge() const { return charge_; }
  inline int quality() const { return quality_; }
  inline int hits_miss_inner() const { return hits_miss_inner_; }


  inline void set_momentum(ThreeVector const& momentum) {
    momentum_ = momentum;
  }
  inline void setPt(double const& pt) { momentum_.SetRho(pt); }
  inline void setEta(double const& eta) { momentum_.SetEta(eta); }
  inline void setPhi(double const& phi) { momentum_.SetPhi(phi); }
  inline void setVx(double const& x) { ref_point_.SetX(x); }
  inline void setVy(double const& y) { ref_point_.SetY(y); }
  inline void setVz(double const& z) { ref_point_.SetZ(z); }
  inline void setHits(int const& hits) {
    hits_ = hits;
  }
  inline void setPixel_hits(int const& pixel_hits) {
    pixel_hits_ = pixel_hits;
  }
  inline void setCharge(int const& charge) { charge_ = charge; }
  inline void setQuality(int const& quality) { quality_ = quality; }
  inline void setHits_miss_inner(int const& hits_miss_inner) {
    hits_miss_inner_ = hits_miss_inner;
  }

 private:
  ThreeVector momentum_;
  Point ref_point_;
  int charge_;
  int hits_;
  int pixel_hits_;
  int quality_;
  int hits_miss_inner_;
};

}
#endif
