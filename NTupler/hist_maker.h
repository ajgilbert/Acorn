#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <vector>
#include <string>
#include <iostream>
#include <math.h>

#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/EventInfo.h"
#include "functions.h"


double hist_maker(std::vector<std::string> input_files, TString outfile_, double cs, TString Dilepton, std::string data_files, int cuts_lev);

double hist_maker_reco(std::vector<std::string> input_files, TString outfile_, double cs, TString Dilepton, std::string data_files, int cuts_lev);
