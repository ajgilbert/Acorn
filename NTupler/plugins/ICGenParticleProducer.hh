#ifndef UserCode_ICHiggsTauTau_ICGenParticleProducer_h
#define UserCode_ICHiggsTauTau_ICGenParticleProducer_h

#include <memory>
#include <vector>
#include <string>
// #include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Rtypes.h"
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/libminifloat.h"
#include "Acorn/NTupler/plugins/ICEventProducer.hh"


// #include "UserCode/ICHiggsTauTau/interface/GenParticle.hh"


class ICGenParticleProducer : public edm::stream::EDProducer<> {
 public:
  explicit ICGenParticleProducer(const edm::ParameterSet &);
  ~ICGenParticleProducer();

  void beginStream(edm::StreamID id) {
    streamid_ = id.value();
  }

  void beginRun(edm::Run const&, edm::EventSetup const&) {
    if (!attachedBranch_) {
      ICEventProducer::getStreamTree(streamid_)->Branch(branch_.c_str(), &particles_);
    }
  }

  template <typename T>
  T reduceMantissaToNbits(const T &f, int bits);

 private:
  // virtual void beginJob();
  virtual void produce(edm::Event &, const edm::EventSetup &);

  std::vector<ROOT::Math::PtEtaPhiMVector> *particles_;
  edm::InputTag input_;
  std::string branch_;
  unsigned streamid_;
  bool attachedBranch_ = false;
  // boost::hash<reco::GenParticle const *> particle_hasher_;

  bool store_mothers_;
  bool store_daughters_;
  bool store_statusFlags_;

  // "drop *"
  // "keep p4=12 pdgId"
  // std::
  struct VarAction {
    bool zeroed = false;
    bool truncate = false;
    int round_to = 0;
  };

  std::map<std::string, VarAction> varActions_;

  template <typename T>
  T SetVar(std::string const& label, T const& arg) {
    auto const& it = varActions_.find(label);
    if (it != varActions_.end()) {
      return it->second.zeroed ? T() : (it->second.truncate ? reduceMantissaToNbits<T>(arg, it->second.round_to) : arg);
    }
    return arg;
  }
  // std::map<std::string, std::function<
};

template <>
double ICGenParticleProducer::reduceMantissaToNbits(const double &f, int bits) {
  uint64_t mask = (0xFFFFFFFFFFFFFFFF >> (52-bits)) << (52-bits);
  union { double flt; uint64_t i64; } conv;
  conv.flt=f;
  conv.i64&=mask;
  return conv.flt;
}

template <>
float ICGenParticleProducer::reduceMantissaToNbits(const float &f, int bits)
{
  uint32_t mask = (0xFFFFFFFF >> (23-bits)) << (23-bits);
  union { float flt; uint32_t i32; } conv;
  conv.flt=f;
  conv.i32&=mask;
  return conv.flt;
}

template <>
ROOT::Math::PtEtaPhiMVector ICGenParticleProducer::reduceMantissaToNbits(const ROOT::Math::PtEtaPhiMVector &v, int bits) {
  if (bits > 0) {
  return ROOT::Math::PtEtaPhiMVector(
      reduceMantissaToNbits<double>(v.pt(), bits),
      reduceMantissaToNbits<double>(v.eta(), bits),
      reduceMantissaToNbits<double>(v.phi(), bits),
      reduceMantissaToNbits<double>(v.M(), bits)
    );
  } else if (bits == -1) {
    union { double x; int16_t y[4]; } conv;
    conv.y[0]  =  MiniFloatConverter::float32to16(float(v.pt()));
    conv.y[1]  =  int16_t(std::round(v.Eta()/6.0f*std::numeric_limits<int16_t>::max()));
    conv.y[2]  =  int16_t(std::round(v.Phi()/3.2f*std::numeric_limits<int16_t>::max()));
    conv.y[3]    =  MiniFloatConverter::float32to16(v.M());
    return ROOT::Math::PtEtaPhiMVector(conv.x, 0., 0., 0.);
  }
  else {
    return v;
  }
}


#endif
