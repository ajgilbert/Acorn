#ifndef Acorn_Analysis_LumiMask_h
#define Acorn_Analysis_LumiMask_h

#include "Acorn/Analysis/interface/TreeEvent.h"
#include "Acorn/Analysis/interface/ModuleBase.h"
#include "Acorn/NTupler/interface/json.hpp"
#include "Acorn/NTupler/interface/RunLumiMap.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include <string>

namespace ac {

  void FillRunLumiMapFromJson(RunLumiMap & rlmap, nlohmann::json const& js);
  nlohmann::json JsonFromRunLumiMap(RunLumiMap const& rlmap);


/**
 * Filters events using a standard CMS luminosity json file
 *
 */
class LumiMask : public ModuleBase {
 private:
  RunLumiMap input_json_;
  RunLumiMap * accept_json_;
  RunLumiMap * reject_json_;
  RunLumiMap * all_json_;
  CLASS_MEMBER(LumiMask, std::string, input_file)
  CLASS_MEMBER(LumiMask, fwlite::TFileService *, fs)

 public:
  LumiMask(std::string const& name);
  virtual ~LumiMask();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent* event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
};
}

#endif
