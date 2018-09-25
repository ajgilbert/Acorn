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
  inline Vector vector() const {
    return Vector(ROOT::Math::PtEtaPhiMVector(pt(), eta(), phi(), 0.13957018));
  }

  inline Point const& ref_point() const { return ref_point_; }

  inline std::size_t id() const { return id_; }
  inline double pt() const { return momentum_.Rho(); }
  inline double energy() const {return momentum_.r(); }
  inline double eta() const { return momentum_.Eta(); }
  inline double phi() const { return momentum_.Phi(); }
  inline double vx() const { return ref_point_.x(); } 
  inline double vy() const { return ref_point_.y(); }
  inline double vz() const { return ref_point_.z(); }
  inline double normalized_chi2() const { return normalized_chi2_; }
  inline int hits() const { return pixel_hits_;}
  //inline double dxy(Point const& point) const {}
  //inline double dz(Point const& point) const {}
  inline int charge() const { return charge_; }
  inline int16_t algorithm() const { return algorithm_; }
  inline double pt_err() const { return pt_err_; }
  inline int quality() const { return quality_; }
  inline int hits_miss_inner() const { return hits_miss_inner_; }


  inline void set_momentum(ThreeVector const& momentum) {
    momentum_ = momentum;
  }
  inline void setId(std::size_t const& id) { id_ = id;}
  inline void setPt(double const& pt) { momentum_.SetRho(pt); }
  inline void setEta(double const& eta) { momentum_.SetEta(eta); }
  inline void setPhi(double const& phi) { momentum_.SetPhi(phi); }
  inline void setVx(double const& x) { ref_point_.SetX(x); }
  inline void setVy(double const& y) { ref_point_.SetY(y); }
  inline void setVz(double const& z) { ref_point_.SetZ(z); }
  inline void setNormalized_chi2(double const& normalized_chi2) {
    normalized_chi2_ = normalized_chi2;
  }
  inline void setHits(int const& hits) {
    hits_ = hits;
  }
  inline void setPixel_hits(int const& pixel_hits) {
    pixel_hits_ = pixel_hits;
  }
  inline void setCharge(int const& charge) { charge_ = charge; }
  inline void setAlgorithm(int16_t const& algorithm) {
    algorithm_ = algorithm;
  }
  inline void setPt_err(double const& pt_err) { pt_err_ = pt_err; }
  inline void setQuality(int const& quality) { quality_ = quality; }
  inline void setHits_miss_inner(int const& hits_miss_inner) {
    hits_miss_inner_ = hits_miss_inner;
  }

 private:
  ThreeVector momentum_;
  Point ref_point_;
  std::size_t id_;
  int charge_;
  double normalized_chi2_;
  int hits_;
  int pixel_hits_;
  int16_t algorithm_;
  double pt_err_;
  int quality_;
  int hits_miss_inner_;
};

}
#endif
