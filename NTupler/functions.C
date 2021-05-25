#include "functions.h"

bool isBJet(ac::Candidate const& jet, ac::GenParticle const& bquark){
        double delta_R_max = .3;
        TLorentzVector jet_lv, b_lv;
        jet_lv.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.M());
        b_lv.SetPtEtaPhiM(bquark.pt(), bquark.eta(), bquark.phi(), bquark.M());
        bool retval = abs(jet_lv.DeltaR(b_lv))< delta_R_max;
        //std::cout<<"isBjet = " << retval << std::endl;
        return retval;
}


//functions to check overlap.  If overlap -> return true
double check_overlap_part(ac::GenParticle const& part_1, ac::GenParticle const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_jet(ac::Candidate const& jet_1, ac::GenParticle const* part_2){
        TLorentzVector jet_lv_1, part_lv_2;
        jet_lv_1.SetPtEtaPhiM(jet_1.pt(), jet_1.eta(), jet_1.phi(), jet_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((jet_lv_1.Phi() - part_lv_2.Phi())*(jet_lv_1.Phi() - part_lv_2.Phi()) + (jet_lv_1.Eta() - part_lv_2.Eta())*(jet_lv_1.Eta() - part_lv_2.Eta())));// < delta_R_max;
        return overlap;
}

double check_overlap_bjet(ac::Candidate const& jet_1, ac::Candidate const* jet_2){
        TLorentzVector jet_lv_1, jet_lv_2;
        jet_lv_1.SetPtEtaPhiM(jet_1.pt(), jet_1.eta(), jet_1.phi(), jet_1.M());
        jet_lv_2.SetPtEtaPhiM(jet_2->pt(), jet_2->eta(), jet_2->phi(), jet_2->M());
        double overlap = abs(sqrt((jet_lv_1.Phi() - jet_lv_2.Phi())*(jet_lv_1.Phi() - jet_lv_2.Phi()) + (jet_lv_1.Eta() - jet_lv_2.Eta())*(jet_lv_1.Eta() - jet_lv_2.Eta())));// < delta_R_max;
        return overlap;
}

double check_overlap(ac::Candidate const& jet_1, ac::GenParticle const& part_2){
        TLorentzVector jet_lv_1, part_lv_2;
        jet_lv_1.SetPtEtaPhiM(jet_1.pt(), jet_1.eta(), jet_1.phi(), jet_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());
        double overlap = abs(sqrt((jet_lv_1.Phi() - part_lv_2.Phi())*(jet_lv_1.Phi() - part_lv_2.Phi()) + (jet_lv_1.Eta() - part_lv_2.Eta())*(jet_lv_1.Eta() - part_lv_2.Eta())));// < delta_R_max;
        return overlap;
}

double check_overlap_pe(ac::Photon const& part_1, ac::Electron const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_pm(ac::Photon const& part_1, ac::Muon const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_bp(ac::PFJet const& part_1, ac::Photon const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_be(ac::PFJet const& part_1, ac::Electron const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_bm(ac::PFJet const& part_1, ac::Muon const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_jp(ac::PFJet const& part_1, ac::Photon const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_je(ac::PFJet const& part_1, ac::Electron const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_jm(ac::PFJet const& part_1, ac::Muon const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}

double check_overlap_jb(ac::PFJet const& part_1, ac::PFJet const* part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2->pt(), part_2->eta(), part_2->phi(), part_2->M());
        double overlap = abs(sqrt((part_lv_1.Phi() - part_lv_2.Phi())*(part_lv_1.Phi() - part_lv_2.Phi()) + (part_lv_1.Eta() - part_lv_2.Eta())*(part_lv_1.Eta() - part_lv_2.Eta()))); // < delta_R_max;
        return overlap;
}
















//Functions used for sorting by descending pT
bool compare_pt_parts(ac::GenParticle const& part_1, ac::GenParticle const& part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

        return (part_lv_1.Pt()> part_lv_2.Pt());
}

bool compare_pt_jets(ac::Candidate const& part_1, ac::Candidate const& part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

        return (part_lv_1.Pt()> part_lv_2.Pt());
}

bool compare_pt_e(ac::Electron const& part_1, ac::Electron const& part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

        return (part_lv_1.Pt()> part_lv_2.Pt());
}

bool compare_pt_m(ac::Muon const& part_1, ac::Muon const& part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

        return (part_lv_1.Pt()> part_lv_2.Pt());
}

bool compare_pt_p(ac::Photon const& part_1, ac::Photon const& part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

        return (part_lv_1.Pt()> part_lv_2.Pt());
}

bool compare_pt_j(ac::PFJet const& part_1, ac::PFJet const& part_2){
        TLorentzVector part_lv_1, part_lv_2;
        part_lv_1.SetPtEtaPhiM(part_1.pt(), part_1.eta(), part_1.phi(), part_1.M());
        part_lv_2.SetPtEtaPhiM(part_2.pt(), part_2.eta(), part_2.phi(), part_2.M());

        return (part_lv_1.Pt()> part_lv_2.Pt());
}





