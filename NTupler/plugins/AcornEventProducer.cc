#include "Acorn/NTupler/plugins/AcornEventProducer.h"
#include <memory>
#include <iostream>
#include <cstdio>
#include "TTree.h"
#include "Compression.h"
#include "TFileMerger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

std::map<unsigned int, TTree *> AcornEventProducer::treeMap_;

AcornEventProducer::AcornEventProducer(const edm::ParameterSet& /*config*/) { }

AcornEventProducer::~AcornEventProducer() {}

void AcornEventProducer::produce(edm::StreamID id, edm::Event&, const edm::EventSetup&) const {
  auto streamProd = this->streamCache(id);
  streamProd->t->Fill();
}

void AcornEventProducer::beginJob() {}

void AcornEventProducer::endJob() {
  std::cout << "AcornEventProducer - merging files...";
  TFileMerger merger(false);
  merger.OutputFile("EventTree.root", true, ROOT::CompressionSettings(ROOT::kLZMA, 5));
  for (auto const& infile : outputFiles_) {
    merger.AddFile(infile);
  }
  merger.Merge();
  std::cout << " done" << std::endl;
  for (auto const& infile : outputFiles_) {
    std::remove(infile.Data());
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(AcornEventProducer);
