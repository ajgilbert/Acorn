#include "Acorn/NTupler/plugins/ICEventProducer.hh"
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// #include "UserCode/ICHiggsTauTau/interface/StaticTree.hh"


std::map<unsigned int, TTree *> ICEventProducer::treeMap_;


ICEventProducer::ICEventProducer(const edm::ParameterSet& /*config*/) {
  // edm::Service<TFileService> fs;
  // ic::StaticTree::tree_ = fs->make<TTree>("EventTree", "EventTree");
}

ICEventProducer::~ICEventProducer() {}

void ICEventProducer::produce(edm::StreamID id, edm::Event&, const edm::EventSetup&) const {
  auto streamProd = this->streamCache(id);
  streamProd->t->Fill();
  // ic::StaticTree::tree_->Fill();
  // ++processed_;
  // if (processed_ == 500) ic::StaticTree::tree_->OptimizeBaskets();
}

void ICEventProducer::beginJob() {}

void ICEventProducer::endJob() {
  std::cout << "ICEventProducer - merging files...";
  TFileMerger merger(false);
  merger.OutputFile("EventTree.root", true, ROOT::CompressionSettings(ROOT::kLZMA, 5));
  for (auto const& infile : outputFiles_) {
    merger.AddFile(infile);
  }
  merger.Merge();
  std::cout << " done" << std::endl;
  for (auto const& infile : outputFiles_) {
    std::remove(infile.Data());
    // std::remove(infile.)
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(ICEventProducer);
