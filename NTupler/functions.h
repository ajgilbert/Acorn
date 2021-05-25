#include "Acorn/NTupler/interface/Candidate.h"
#include "Acorn/NTupler/interface/GenParticle.h"
#include "Acorn/NTupler/interface/Met.h"
#include "Acorn/NTupler/interface/EventInfo.h"

#include "Acorn/NTupler/interface/Electron.h"
#include "Acorn/NTupler/interface/Muon.h"
#include "Acorn/NTupler/interface/Photon.h"
#include "Acorn/NTupler/interface/PFJet.h"


#include <math.h>

bool isBJet(ac::Candidate const& jet, ac::GenParticle const& bquark);
//functions to check overlap.  If overlap -> return true

double check_overlap_part(ac::GenParticle const& part_1, ac::GenParticle const* part_2);

double check_overlap_jet(ac::Candidate const& jet_1, ac::GenParticle const* part_2);

double check_overlap_bjet(ac::Candidate const& jet_1, ac::Candidate const* jet_2);

double check_overlap(ac::Candidate const& jet_1, ac::GenParticle const& part_2);

double check_overlap(ac::Candidate const& jet_1, ac::GenParticle const& part_2);

double check_overlap_pe(ac::Photon const& part_1, ac::Electron const* part_2);

double check_overlap_pm(ac::Photon const& part_1, ac::Muon const* part_2);

double check_overlap_bp(ac::PFJet const& part_1, ac::Photon const* part_2);

double check_overlap_be(ac::PFJet const& part_1, ac::Electron const* part_2);

double check_overlap_bm(ac::PFJet const& part_1, ac::Muon const* part_2);

double check_overlap_jp(ac::PFJet const& part_1, ac::Photon const* part_2);

double check_overlap_je(ac::PFJet const& part_1, ac::Electron const* part_2);

double check_overlap_jm(ac::PFJet const& part_1, ac::Muon const* part_2);

double check_overlap_jb(ac::PFJet const& part_1, ac::PFJet const* part_2);






//Functions used for sorting by descending pT
bool compare_pt_parts(ac::GenParticle const& part_1, ac::GenParticle const& part_2);

bool compare_pt_jets(ac::Candidate const& part_1, ac::Candidate const& part_2);

bool compare_pt_e(ac::Electron const& part_1, ac::Electron const& part_2);

bool compare_pt_m(ac::Muon const& part_1, ac::Muon const& part_2);

bool compare_pt_p(ac::Photon const& part_1, ac::Photon const& part_2);

bool compare_pt_j(ac::PFJet const& part_1, ac::PFJet const& part_2);


