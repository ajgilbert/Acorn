#ifndef UserCode_ICHiggsTauTau_RequestByDeltaR_h
#define UserCode_ICHiggsTauTau_RequestByDeltaR_h

#include <memory>
#include <vector>
#include "boost/functional/hash.hpp"
#include "boost/format.hpp"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

template <class T, class U>
class RequestByDeltaR : public edm::EDProducer {
 public:
  explicit RequestByDeltaR(const edm::ParameterSet &);
  ~RequestByDeltaR();

 private:
  virtual void beginJob();
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endJob();

  edm::EDGetTokenT<edm::View<T>> input_;
  edm::EDGetTokenT<edm::View<U>> reference_;
  double dr_;

  typedef std::vector<T> Vec;
  typedef edm::RefVector<Vec> RefVectorVec;
  typedef edm::Ref<Vec> RefVec;
};

// =============================
// Template class implementation
// =============================
template <class T, class U>
RequestByDeltaR<T,U>::RequestByDeltaR(const edm::ParameterSet& config)
    : input_(consumes<edm::View<T>>(config.getParameter<edm::InputTag>("src"))),
      reference_(consumes<edm::View<U>>(config.getParameter<edm::InputTag>("reference"))),
      dr_(config.getParameter<double>("deltaR")) {
  produces<RefVectorVec>();
}

template <class T, class U>
RequestByDeltaR<T, U>::~RequestByDeltaR() { }

// =============
// Main producer
// =============
template <class T, class U>
void RequestByDeltaR<T, U>::produce(edm::Event& event,
                                 const edm::EventSetup& setup) {

  edm::Handle<edm::View<T> > ref_handle;
  event.getByToken(reference_, ref_handle);
  edm::Handle<edm::View<U> > in_handle;
  event.getByToken(input_, in_handle);

  std::unique_ptr<RefVectorVec> product(new RefVectorVec());

  for (unsigned i = 0; i < in_handle->size(); ++i) {
    for (unsigned j = 0; j < ref_handle->size(); ++j) {
      double dr = reco::deltaR(in_handle->at(i), ref_handle->at(j));
      if (dr <= dr_) {
        product->push_back(in_handle->refAt(i).template castTo<RefVec>());
        break;
      }
    }
  }
  event.put(std::move(product));
}

template <class T, class U>
void RequestByDeltaR<T, U>::beginJob() {}

template <class T, class U>
void RequestByDeltaR<T, U>::endJob() {}

#endif
