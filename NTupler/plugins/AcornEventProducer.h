#ifndef Acorn_NTupler_EventProducer_h
#define Acorn_NTupler_EventProducer_h

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

struct StreamProducer {
  TFile *f = nullptr;
  TTree *t = nullptr;
};


class AcornEventProducer : public edm::global::EDProducer<edm::StreamCache<StreamProducer>> {
 public:
  explicit AcornEventProducer(const edm::ParameterSet&);
  ~AcornEventProducer();

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
    return treeMap_.at(id);
  }

  template <class T>
  static void AddBranch(unsigned int id, std::string name, T* ptr) {
    auto treeptr = getStreamTree(id);
    m_queue.pushAndWait([&]() {
      treeptr->Branch(name.c_str(), &ptr);
    });
  }

 private:
  virtual void beginJob();
  virtual void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const;
  virtual void endJob();

  static std::map<unsigned int, TTree *> treeMap_;
  mutable std::vector<TString> outputFiles_;

  static edm::SerialTaskQueue m_queue;
};

#endif
