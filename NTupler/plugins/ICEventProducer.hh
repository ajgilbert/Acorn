#ifndef UserCode_ICHiggsTauTau_ICEventProducer_h
#define UserCode_ICHiggsTauTau_ICEventProducer_h

#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "Compression.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Concurrency/interface/SerialTaskQueue.h"

/**
 * @brief Handles the creation of the ntuple output, and must always be included
 *after the other IC object producers.
 *
 * **Example usage**
 * @snippet python/default_producers_cfi.py Event
 */


struct StreamProducer {
  TFile *f = nullptr;
  TTree *t = nullptr;
};


class ICEventProducer : public edm::global::EDProducer<edm::StreamCache<StreamProducer>> {
 public:
  explicit ICEventProducer(const edm::ParameterSet&);
  ~ICEventProducer();

  std::unique_ptr<StreamProducer> beginStream(edm::StreamID id) const {
      auto streamProd = std::make_unique<StreamProducer>();
      m_queue.pushAndWait([&]() {
        TString filename = TString::Format("EventTree_%i.root", id.value());
        streamProd->f = new TFile(filename, "RECREATE");
        streamProd->f->SetCompressionSettings(ROOT::CompressionSettings(ROOT::kLZMA, 5));
        outputFiles_.push_back(filename);
        streamProd->t = new TTree("EventTree", "EventTree");
        this->treeMap_[id.value()] = streamProd->t;
      });
      return streamProd;
  }

  void endStream(edm::StreamID id) const {
      auto streamProd = this->streamCache(id);
      m_queue.pushAndWait([&]() {
        streamProd->f->cd();
        streamProd->t->Write();
        streamProd->f->Close();
      });
  }

  static TTree* getStreamTree(unsigned int id) {
    return treeMap_[id];
  }

 private:
  virtual void beginJob();
  virtual void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const;
  virtual void endJob();

  static std::map<unsigned int, TTree *> treeMap_;
  mutable std::vector<TString> outputFiles_;

  mutable edm::SerialTaskQueue m_queue;
};

#endif
