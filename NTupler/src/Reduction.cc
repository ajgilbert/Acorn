#include <memory>
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "DataFormats/PatCandidates/interface/libminifloat.h"
#include "Acorn/NTupler/interface/Reduction.h"

template <>
double reduceMantissaToNbitsRounding(const double &f, int bits) {
  uint32_t shift = 23 - bits;
  uint32_t mask = (0xFFFFFFFF >> (shift)) << (shift);
  uint32_t test = 1 << (shift - 1);
  uint32_t maxn = (1 << bits) - 2;
  constexpr uint32_t low23 = (0x007FFFFF);  // mask to keep lowest 23 bits = mantissa
  constexpr uint32_t hi9 = (0xFF800000);    // mask to keep highest 9 bits = the rest
  union {
    float flt;
    uint32_t i32;
  } conv;
  conv.flt = f;
  if (conv.i32 & test) {  // need to round
    uint32_t mantissa = (conv.i32 & low23) >> shift;
    if (mantissa < maxn) mantissa++;
    conv.i32 = (conv.i32 & hi9) | (mantissa << shift);
  } else {
    conv.i32 &= mask;
  }
  return conv.flt;
}

template <>
float reduceMantissaToNbitsRounding(const float &f, int bits) {
  uint64_t shift = 52 - bits;
  uint64_t mask = (0xFFFFFFFFFFFFFFFF >> (shift)) << (shift);
  uint64_t test = 1 << (shift - 1);
  uint64_t maxn = (1 << bits) - 2;
  constexpr uint64_t low52 = (0x000FFFFFFFFFFFFF);  // mask to keep lowest 23 bits = mantissa
  constexpr uint64_t hi12 = (0xFFF0000000000000);   // mask to keep highest 9 bits = the rest
  union {
    double flt;
    uint64_t i64;
  } conv;
  conv.flt = f;
  if (conv.i64 & test) {  // need to round
    uint64_t backup = conv.i64;
    conv.i64 &= mask;
    conv.i64 = backup;
    uint64_t mantissa = (conv.i64 & low52) >> shift;
    if (mantissa < maxn) mantissa++;
    conv.i64 = (conv.i64 & hi12) | (mantissa << shift);

  } else {
    conv.i64 &= mask;
  }
  return conv.flt;
}

template <>
double reduceMantissaToNbits(const double &f, int bits) {
  uint64_t mask = (0xFFFFFFFFFFFFFFFF >> (52 - bits)) << (52 - bits);
  union {
    double flt;
    uint64_t i64;
  } conv;
  conv.flt = f;
  conv.i64 &= mask;
  return conv.flt;
}

template <>
float reduceMantissaToNbits(const float &f, int bits) {
  uint32_t mask = (0xFFFFFFFF >> (23 - bits)) << (23 - bits);
  union {
    float flt;
    uint32_t i32;
  } conv;
  conv.flt = f;
  conv.i32 &= mask;
  return conv.flt;
}

template <>
double Reduce(const double &f, int bits) {
  return reduceMantissaToNbitsRounding<double>(f, bits);
}

template <>
float Reduce(const float &f, int bits) {
  return reduceMantissaToNbitsRounding<float>(f, bits);
}

template <>
ROOT::Math::PtEtaPhiMVector Reduce(const ROOT::Math::PtEtaPhiMVector &v, int bits) {
  if (bits > 0) {
    return ROOT::Math::PtEtaPhiMVector(reduceMantissaToNbitsRounding<double>(v.pt(), bits),
                                       reduceMantissaToNbitsRounding<double>(v.eta(), bits),
                                       reduceMantissaToNbitsRounding<double>(v.phi(), bits),
                                       reduceMantissaToNbitsRounding<double>(v.M(), bits));
  } else if (bits == -1) {  // Instead store the P4 just like in a PackedCandidate
    union {
      double x;
      int16_t y[4];
    } conv;
    conv.y[0] = MiniFloatConverter::float32to16(float(v.pt()));
    conv.y[1] = int16_t(std::round(v.Eta() / 6.0f * std::numeric_limits<int16_t>::max()));
    conv.y[2] = int16_t(std::round(v.Phi() / 3.2f * std::numeric_limits<int16_t>::max()));
    conv.y[3] = MiniFloatConverter::float32to16(v.M());

    // All four 16-bit floats are now stored in a single double
    return ROOT::Math::PtEtaPhiMVector(conv.x, 0., 0., 0.);
  } else {
    return v;
  }
}

template <>
ROOT::Math::XYZPoint Reduce(const ROOT::Math::XYZPoint &v, int bits) {
  if (bits > 0) {
    return ROOT::Math::XYZPoint(reduceMantissaToNbitsRounding<double>(v.x(), bits),
                                reduceMantissaToNbitsRounding<double>(v.y(), bits),
                                reduceMantissaToNbitsRounding<double>(v.z(), bits));
  } else {
    return v;
  }
}

