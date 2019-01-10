#ifndef Acorn_NTupler_EventInfo_hh
#define Acorn_NTupler_EventInfo_hh

#include <iostream>
#include <map>
#include <string>
#include <bitset>
// #include "Acorn/NTupler/interface/city.h"
#include "Rtypes.h"

namespace ac {

class EventInfo {
 private:
 public:
  EventInfo();
  virtual ~EventInfo();
  virtual void Print(unsigned detail=0) const;

  /// If event is real data returns `true`, otherwise `false`
  inline bool isRealData() const { return isRealData_; }

  /// Event number
  inline unsigned long long event() const { return event_; }

  /// Run number
  inline int run() const { return run_; }

  /// Lumisection number
  inline int luminosityBlock() const { return luminosityBlock_; }

  /// Bunch crossing number
  inline int bunchCrossing() const { return bunchCrossing_; }

  inline double nominalGenWeight() const { return nominalGenWeight_; }

  inline std::vector<double> const& genWeights() const { return genWeights_; }

  inline double nominalLHEWeight() const { return nominalLHEWeight_; }

  inline std::map<unsigned, double> const& lheWeights() const { return lheWeights_; }

  inline std::bitset<32> metfilters() const { return std::bitset<32>(metfilters_); }

  inline std::vector<double> const& userDoubles() const { return userDoubles_; }

  // /// Energy density used for the jet energy corrections in this event
  // inline double jet_rho() const { return jet_rho_; }

  // /// Energy density used for calculating lepton isolation in this event
  // inline double lepton_rho() const { return lepton_rho_; }

  // /// Generator level HT, used for combining HT binned samples with inclusive samples
  // inline double gen_ht() const { return gen_ht_; }

  // /// Generator level M_ll
  // inline double gen_mll() const { return gen_mll_; }

  // /// Number of outgoing partons at generator level, used for combining n-jet binned samples with
  // inclusive samples inline unsigned n_outgoing_partons() const { return n_outgoing_partons_; }

  /// Number of reconstructed vertices passing some baseline quality
  /// requirements
  // inline unsigned good_vertices() const { return good_vertices_; }
  /**@}*/

  inline void setIsRealData(bool const& isRealData) { isRealData_ = isRealData; }

  inline void setEvent(unsigned long long const& event) { event_ = event; }

  inline void setRun(int const& run) { run_ = run; }

  inline void setLuminosityBlock(int const& luminosityBlock) { luminosityBlock_ = luminosityBlock; }

  inline void setBunchCrossing(int const& bunchCrossing) { bunchCrossing_ = bunchCrossing; }

  inline void setNominalGenWeight(double const& nominalGenWeight) {
    nominalGenWeight_ = nominalGenWeight;
  }

  inline void setGenWeights(std::vector<double> const& genWeights) { genWeights_ = genWeights; }

  inline void setNominalLHEWeight(double const& nominalLHEWeight) {
    nominalLHEWeight_ = nominalLHEWeight;
  }

  inline void setLHEWeight(unsigned const& id, double const& weight) { lheWeights_[id] = weight; }

  inline void setMetFilters(std::bitset<32> const& metfilters) { metfilters_ = unsigned(metfilters.to_ulong()); }

  inline void setUserDoubles(std::vector<double> const& userDoubles) { userDoubles_ = userDoubles; }

  inline void setWeight(std::string const& label, double const& wt, bool const& enabled = true) {
    weights_[label] = std::make_pair(enabled, wt);
  }

  inline double totalWeight() const {
    double total = 1.0;
    for (auto const& it : weights_) {
      if (it.second.first) {
        total *= it.second.second;
      }
    }
    return total;
  }

  // /// @copybrief jet_rho()
  // inline void set_jet_rho(double const& jet_rho) { jet_rho_ = jet_rho; }

  // /// @copybrief lepton_rho()
  // inline void set_lepton_rho(double const& lepton_rho) {
  //   lepton_rho_ = lepton_rho;
  // }

  // /// @copybrief gen_ht()
  // inline void set_gen_ht(double const& gen_ht) { gen_ht_ = gen_ht; }

  // /// @copybrief gen_mll()
  // inline void set_gen_mll(double const& gen_mll) { gen_mll_ = gen_mll; }

  // /// @copybrief n_outgoing_partons()
  // inline void set_n_outgoing_partons(unsigned const& n_outgoing_partons) { n_outgoing_partons_ =
  // n_outgoing_partons; }

  // /// @copybrief good_vertices()
  // inline void set_good_vertices(unsigned const& good_vertices) {
  //   good_vertices_ = good_vertices;
  // }
  /**@}*/

 private:
  // General event properties
  bool isRealData_;
  unsigned long long event_;
  int run_;
  int luminosityBlock_;
  int bunchCrossing_;

  // Generator weight (should be present for all MC)
  double nominalGenWeight_;
  // Will include PYTHIA weights, e.g. for PS+UE variations. Given relative to nominalGenWeight_
  std::vector<double> genWeights_;

  // LHE weights (present only for samples with external generators)
  double nominalLHEWeight_;
  std::map<unsigned, double> lheWeights_;  // given relative to nominalLHEWeight_

  std::map<std::string, std::pair<bool, double>> weights_;

  unsigned metfilters_;

  std::vector<double> userDoubles_;

  // double jet_rho_;
  // double lepton_rho_;
  // double gen_ht_;
  // unsigned n_outgoing_partons_;
  // double gen_mll_;
  // SDMap weights_;
  // SBMap weight_status_;
  // unsigned good_vertices_;
  // TBMap filters_;
};
}  // namespace ac
#endif
