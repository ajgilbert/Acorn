#include "Acorn/Analysis/interface/StandaloneReweight.h"
#include "boost/lexical_cast.hpp"
#include "Python.h"
#include "boost/python.hpp"
#include <fstream>
#include <vector>
#include <iostream>


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

std::vector<double> StandaloneReweight::TransformWeights(std::vector<double> const& in) const {
  int verb = 0;
  if (verb > 0) std::cout << "-- Have " << in.size() << " weights\n";
  std::cout.precision(9);
  std::vector<double> out(in.size(), 0.);
  // unsigned N = (in.size() - 2) / 2;

  // We subtract 1 here, not 2
  unsigned N = (-3 + int(sqrt(9.0 + 8.0 * (double(in.size()) - 1.0)) + 0.5)) / 2;
  for (unsigned i = 0; i < in.size(); ++i) {
    if (verb > 0) std::cout << " - " << in[i] << "\n";
    out[i] = in[i];
  }
  for (unsigned ip = 0; ip < N; ++ip) {
    double s0 = in[0];
    double s1 = in[ip * 2 + 1];
    double s2 = in[ip * 2 + 2];
    if (verb > 0) std::cout << " -- Doing " << ip << "\n";
    if (verb > 0) std::cout << s0 << "\t" << s1 << "\t" << s2 << "\n";
    s1 -= s0;
    s2 -= s0;
    if (verb > 0) std::cout << " - subtract s0: " << s1 << "\t" << s2 << "\n";
    double Ai = 4. * s1 - s2;
    double Bii = s2 - Ai;
    if (verb > 0) std::cout << " - Result: " << Ai << "\t" << Bii << "\n";
    out[ip * 2 + 1] = Ai;
    out[ip * 2 + 2] = Bii;
  }
  unsigned crossed_offset = 1 + 2 * N;
  unsigned c_counter = 0;
  for (unsigned ix = 0; ix < N; ++ix) {
    for (unsigned iy = ix + 1; iy < N; ++iy) {
      if (verb > 0) std::cout << " -- Doing " << ix << "\t" << iy << "\t[" << (crossed_offset + c_counter) << "]\n";
      double s = in[crossed_offset + c_counter];
      double sm = in[0];
      double sx = out[ix * 2 + 1];
      double sy = out[iy * 2 + 1];
      double sxx = out[ix * 2 + 2];
      double syy = out[iy * 2 + 2];
      s -= (sm + sx + sy + sxx + syy);
      out[crossed_offset + c_counter] = s;
      if (verb > 0) std::cout << " - Result: " << s << "\n";
      ++c_counter;
    }
  }
  // for (unsigned i = 0; i < in.size(); ++i) {
  //   in[i] = out[i];
  // }
  return out;
}