#include "Acorn/Analysis/interface/StandaloneReweight.h"
#include "boost/lexical_cast.hpp"
#include "Python.h"
#include "boost/python.hpp"
#include <fstream>
#include <vector>


StandaloneReweight::StandaloneReweight(std::string const& rw_module, std::string const& rw_pack) {
  namespace py = boost::python;
  Py_Initialize();
  py::object mod = py::import(rw_module.c_str());
  py::object StandaloneReweight = mod.attr("StandaloneReweight");
  rw_ = StandaloneReweight(rw_pack);
}

StandaloneReweight::~StandaloneReweight() { ; }

std::vector<double> StandaloneReweight::ComputeWeights(
    std::vector<std::vector<double>> const& parts, std::vector<int> const& pdgs,
    std::vector<int> const& hels, std::vector<int> const& stats, double alphas, bool doHelicity,
    bool verb) {
  boost::python::list py_parts;
  for (auto const& part : parts) {
    py_parts.append(MakeList(part));
  }

  boost::python::list py_pdgs = MakeList(pdgs);
  boost::python::list py_hels = MakeList(hels);
  boost::python::list py_stats = MakeList(stats);

  boost::python::list py_res(
      rw_.attr("ComputeWeights")(py_parts, py_pdgs, py_hels, py_stats, alphas, doHelicity, verb));
  std::vector<double> res;
  for (unsigned i = 0; i < boost::python::len(py_res); ++i) {
    res.push_back(boost::python::extract<double>(py_res[i]));
  }
  return res;
}
