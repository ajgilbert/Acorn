#ifndef Acorn_RunLumiMap_hh
#define Acorn_RunLumiMap_hh
#include <map>
#include <set>
#include <vector>
#include "Rtypes.h"
#include "TNamed.h"
#include "TCollection.h"

class RunLumiMap : public TNamed {
 public:
  RunLumiMap();
  RunLumiMap(const char *name, const char *title);
  virtual ~RunLumiMap();

  Long64_t Merge(TCollection* coll);

  void Add(unsigned run, unsigned ls);
  void Remove(unsigned run, unsigned ls);
  bool InMap(unsigned run, unsigned ls) const;

  std::map<unsigned, std::set<unsigned>> const& GetMap() const;
  void ClearMap();
  void ResetMap(std::map<unsigned, std::set<unsigned>> const& newmap);
  void AppendMap(std::map<unsigned, std::set<unsigned>> const& newmap);
  void PrintMap() const;

  void AppendFromJsonString(std::string const& js_str);
  std::string AsJsonString() const;


  ClassDef(RunLumiMap,1)

 private:
  std::map<unsigned, std::set<unsigned>> lsmap_;
};

#endif
