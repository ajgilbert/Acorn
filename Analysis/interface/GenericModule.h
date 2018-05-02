#ifndef ICHiggsTauTau_Module_GenericModule_h
#define ICHiggsTauTau_Module_GenericModule_h

#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "boost/function.hpp"

#include <string>

namespace ac {

class GenericModule : public ModuleBase {
 private:
  CLASS_MEMBER(GenericModule, boost::function<int(ac::TreeEvent *)>, function)

 public:
  GenericModule(std::string const& name);
  virtual ~GenericModule();

  virtual int PreAnalysis();
  virtual int Execute(ac::TreeEvent* evt);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};

GenericModule::GenericModule(std::string const& name) : ModuleBase(name) {
}

GenericModule::~GenericModule() {
  ;
}

int GenericModule::PreAnalysis() {
  return 0;
}


int GenericModule::Execute(ac::TreeEvent* event) {
  return function_(event);
}


int GenericModule::PostAnalysis() {
  return 0;
}

void GenericModule::PrintInfo() {
  ;
}





}

#endif
