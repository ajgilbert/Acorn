#ifndef Acorn_NTupler_Reduction_h
#define Acorn_NTupler_Reduction_h

#include <memory>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"
#include "DataFormats/PatCandidates/interface/libminifloat.h"

template <typename T>
T reduceMantissaToNbitsRounding(const T &f, int bits);

template <>
double reduceMantissaToNbitsRounding(const double &f, int bits);

template <>
float reduceMantissaToNbitsRounding(const float &f, int bits);

template <typename T>
T reduceMantissaToNbits(const T &f, int bits);

template <>
double reduceMantissaToNbits(const double &f, int bits);

template <>
float reduceMantissaToNbits(const float &f, int bits);

// The default implementation of Reduce just returns the same value
// Below this we specialise for different types
template <typename T>
T Reduce(const T &f, int /*bits*/) {
  return f;
}

template <>
double Reduce(const double &f, int bits);

template <>
float Reduce(const float &f, int bits);

template <>
ROOT::Math::PtEtaPhiMVector Reduce(const ROOT::Math::PtEtaPhiMVector &v, int bits);

template <>
ROOT::Math::XYZPoint Reduce(const ROOT::Math::XYZPoint &v, int bits);

#endif
