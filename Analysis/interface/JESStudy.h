#ifndef Acorn_Analysis_JESStudy_h
#define Acorn_Analysis_JESStudy_h
#include <string>
#include <cstdint>
#include "boost/range/algorithm_ext/erase.hpp"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include "TRandom3.h"
#include "TH2D.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/Analysis/interface/AnalysisTools.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/Analysis/interface/RoccoR.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

namespace ac {

class JESStudy : public ModuleBase {
 private:
  CLASS_MEMBER(JESStudy, fwlite::TFileService*, fs)
  CLASS_MEMBER(JESStudy, unsigned, year)
  CLASS_MEMBER(JESStudy, bool, is_data)

  std::vector<JetCorrectorParameters *> jpars_;
  std::vector<JetCorrectionUncertainty *> juncs_;
  std::vector<std::string> jsrcs_;
  std::vector<TH2D *> jhists_;

 private:
  TH2D *j_nominal;

 public:
  JESStudy(std::string const& name);
  virtual ~JESStudy();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}  // namespace ac

#endif
