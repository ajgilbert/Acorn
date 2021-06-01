#ifndef Acorn_Analysis_StandaloneReweight_h
#define Acorn_Analysis_StandaloneReweight_h

#include <string>
#include <vector>
#include "boost/python.hpp"

class StandaloneReweight {
 private:
  boost::python::object rw_;

  template <class T>
  boost::python::list MakeList(std::vector<T> const& in);

 public:
  StandaloneReweight(std::string const& rw_module, std::string const& rw_pack);
  ~StandaloneReweight();

  std::vector<double> ComputeWeights(std::vector<std::vector<double>> const& parts,
                                     std::vector<int> const& pdgs, std::vector<int> const& hels,
                                     std::vector<int> const& stats, double alphas, bool doHelicity,
                                     bool verb);

  std::vector<double> TransformWeights(std::vector<double> const& in) const;
};

template <class T>
boost::python::list StandaloneReweight::MakeList(std::vector<T> const& in) {
  boost::python::list res;
  for (auto const& arg : in) {
    res.append(arg);
  }
  return res;
}

#endif
