#!/usr/bin/env python
import ROOT
# import json
import argparse
from array import array
from Acorn.Analysis.workspaceTools import *


def GetFromTFile(str):
    f = ROOT.TFile(str.split(':')[0])
    obj = f.Get(str.split(':')[1]).Clone()
    f.Close()
    return obj

# Boilerplate
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.RooWorkspace.imp = getattr(ROOT.RooWorkspace, 'import')
ROOT.TH1.AddDirectory(0)

parser = argparse.ArgumentParser()

parser.add_argument('--era', default='2016',
                    help='The dataset to target')
parser.add_argument('--debug', action='store_true',
                    help='Print debug output')
args = parser.parse_args()

# ROOT.gROOT.LoadMacro("CrystalBallEfficiency.cxx+")

w = ROOT.RooWorkspace('w')

era = args.era


###############################################################################
## Pileup
###############################################################################

if era == '2016':
    mc_edges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
                35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
                67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
                83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]
    mc_probs = [1.78653e-05, 2.56602e-05, 5.27857e-05, 8.88954e-05, 0.000109362,
                0.000140973, 0.000240998, 0.00071209, 0.00130121, 0.00245255,
                0.00502589, 0.00919534, 0.0146697, 0.0204126, 0.0267586, 0.0337697,
                0.0401478, 0.0450159, 0.0490577, 0.0524855, 0.0548159, 0.0559937,
                0.0554468, 0.0537687, 0.0512055, 0.0476713, 0.0435312, 0.0393107,
                0.0349812, 0.0307413, 0.0272425, 0.0237115, 0.0208329, 0.0182459,
                0.0160712, 0.0142498, 0.012804, 0.011571, 0.010547, 0.00959489,
                0.00891718, 0.00829292, 0.0076195, 0.0069806, 0.0062025, 0.00546581,
                0.00484127, 0.00407168, 0.00337681, 0.00269893, 0.00212473, 0.00160208,
                0.00117884, 0.000859662, 0.000569085, 0.000365431, 0.000243565,
                0.00015688, 9.88128e-05, 6.53783e-05, 3.73924e-05, 2.61382e-05,
                2.0307e-05, 1.73032e-05, 1.435e-05, 1.36486e-05, 1.35555e-05,
                1.37491e-05, 1.34255e-05, 1.33987e-05, 1.34061e-05, 1.34211e-05,
                1.34177e-05, 1.32959e-05, 1.33287e-05]
    h_data = GetFromTFile('hvm/inputs/pileup_profile_2016.root:pileup')

if era == '2017':
    mc_edges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
                35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
                67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
                83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]
    mc_probs = [3.39597497605e-05, 6.63688402133e-06, 1.39533611284e-05, 3.64963078209e-05,
                6.00872171664e-05, 9.33932578027e-05, 0.000120591524486, 0.000128694546198,
                0.000361697233219, 0.000361796847553, 0.000702474896113, 0.00133766053707,
                0.00237817050805, 0.00389825605651, 0.00594546732588, 0.00856825906255,
                0.0116627396044, 0.0148793350787, 0.0179897368379, 0.0208723871946,
                0.0232564170641, 0.0249826433945, 0.0262245860346, 0.0272704617569,
                0.0283301107549, 0.0294006137386, 0.0303026836965, 0.0309692426278,
                0.0308818046328, 0.0310566806228, 0.0309692426278, 0.0310566806228,
                0.0310566806228, 0.0310566806228, 0.0307696426944, 0.0300103336052,
                0.0288355370103, 0.0273233309106, 0.0264343533951, 0.0255453758796,
                0.0235877272306, 0.0215627588047, 0.0195825559393, 0.0177296309658,
                0.0160560731931, 0.0146022004183, 0.0134080690078, 0.0129586991411,
                0.0125093292745, 0.0124360740539, 0.0123547104433, 0.0123953922486,
                0.0124360740539, 0.0124360740539, 0.0123547104433, 0.0124360740539,
                0.0123387597772, 0.0122414455005, 0.011705203844, 0.0108187105305,
                0.00963985508986, 0.00827210065136, 0.00683770076341, 0.00545237697118,
                0.00420456901556, 0.00367513566191, 0.00314570230825, 0.0022917978982,
                0.00163221454973, 0.00114065309494, 0.000784838366118, 0.000533204105387,
                0.000358474034915, 0.000238881117601, 0.0001984254989, 0.000157969880198,
                0.00010375646169, 6.77366175538e-05, 4.39850477645e-05, 2.84298066026e-05,
                1.83041729561e-05, 1.17473542058e-05, 7.51982735129e-06, 6.16160108867e-06,
                4.80337482605e-06, 3.06235473369e-06, 1.94863396999e-06, 1.23726800704e-06,
                7.83538083774e-07, 4.94602064224e-07, 3.10989480331e-07, 1.94628487765e-07,
                1.57888581037e-07, 1.2114867431e-07, 7.49518929908e-08, 4.6060444984e-08,
                2.81008884326e-08, 1.70121486128e-08, 1.02159894812e-08]
    h_data = GetFromTFile('hvm/inputs/pileup_profile_2017.root:pileup')

if era == '2018':
    mc_edges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
                39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62,
                63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
                75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,
                87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]
    mc_probs = [4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05,
                3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473,
                0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138,
                0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411,
                0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554,
                0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895,
                0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877,
                0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612,
                0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551,
                0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934,
                0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915,
                0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932,
                0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885,
                0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012,
                0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05,
                2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06,
                3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07,
                5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07,
                1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08,
                6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08]
    h_data = GetFromTFile('hvm/inputs/pileup_profile_2018.root:pileup')

print len(mc_edges), len(mc_probs)

h_mc = ROOT.TH1D('mc_pileup', '', len(mc_edges), mc_edges[0], mc_edges[-1] + 1)
for i in range(len(mc_edges)):
    if i >= len(mc_probs):
        h_mc.SetBinContent(i + 1, 0.)
    else:
        h_mc.SetBinContent(i + 1, mc_probs[i])

# Fix for very large weights in 2018
# if era == '2018':
#     restrict = 65
#     h_mc_new = ROOT.TH1D('mc_pileup', '', restrict, 0., float(restrict))
#     h_data_new = ROOT.TH1D('data_pileup', '', restrict, 0., float(restrict))
#     for i in range(restrict):
#         h_mc_new.SetBinContent(i + 1, h_mc.GetBinContent(i + 1))
#         h_data_new.SetBinContent(i + 1, h_data.GetBinContent(i + 1))
#     h_mc = h_mc_new
#     h_data = h_data_new

h_data.Scale(1. / h_data.Integral())
h_mc.Scale(1. / h_mc.Integral())

SafeWrapHist(w, ['pu_int'], h_data, name='pileup_data')
SafeWrapHist(w, ['pu_int'], h_mc, name='pileup_mc')

w.factory('expr::pileup_ratio("@0/@1", pileup_data, pileup_mc)')

if args.debug:
    print '>> Debug pileup weights:'
    x = 0.5
    for i in xrange(100):
        w.var('pu_int').setVal(x)
        print 'pu_int = %f, pileup_ratio = %f' % (x, w.function('pileup_ratio').getVal())
        x += 1.0

###############################################################################
## Muons
###############################################################################
loc = 'hvm/inputs/muons/%s' % era
if era == '2016':
    histsToWrap = [
        (loc + '/EfficienciesAndSF_RunBtoF.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_bf_ratio'),
        (loc + '/EfficienciesAndSF_RunBtoF.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/pt_abseta_DATA', 'm_trg_bf_data_eff'),
        (loc + '/EfficienciesAndSF_RunBtoF.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/pt_abseta_MC', 'm_trg_mc_eff'),
        (loc + '/EfficienciesAndSF_RunGtoH.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_gh_ratio'),
        (loc + '/EfficienciesAndSF_RunGtoH.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/pt_abseta_DATA', 'm_trg_gh_data_eff')
    ]
    histsToWrapNewStyle = [
        (loc + '/RunBCDEF_SF_ID.root:NUM_MediumID_DEN_genTracks_eta_pt', 'm_id_bf_ratio'),
        (loc + '/RunGH_SF_ID.root:NUM_MediumID_DEN_genTracks_eta_pt', 'm_id_gh_ratio'),
        (loc + '/RunBCDEF_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_eta_pt', 'm_iso_bf_ratio'),
        (loc + '/RunGH_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_eta_pt', 'm_iso_gh_ratio')
    ]

    #muon_trk_eff_hist = TGraphAsymmErrorsToTH1D(GetFromTFile(loc + '/track/Tracking_EfficienciesAndSF_BCDEFGH.root:ratio_eff_eta3_dr030e030_corr'))
    #SafeWrapHist(w, ['m_eta'], muon_trk_eff_hist, name='m_trk_ratio')

if era == '2017':
    histsToWrap = [
        (loc + '/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root:IsoMu27_PtEtaBins/pt_abseta_ratio', 'm_trg_ratio'),
        (loc + '/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root:IsoMu27_PtEtaBins/efficienciesDATA/pt_abseta_DATA', 'm_trg_data_eff'),
        (loc + '/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root:IsoMu27_PtEtaBins/efficienciesMC/pt_abseta_MC', 'm_trg_mc_eff'),
        (loc + '/RunBCDEF_SF_ID.root:NUM_MediumID_DEN_genTracks_pt_abseta', 'm_id_ratio'),
        (loc + '/RunBCDEF_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_pt_abseta', 'm_iso_ratio'),
    ]
    histsToWrapNewStyle = []

if era == '2018':
    histsToWrap = [
        (loc + '/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root:IsoMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_a_ratio'),
        (loc + '/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root:IsoMu24_PtEtaBins/efficienciesDATA/pt_abseta_DATA', 'm_trg_a_data_eff'),
        (loc + '/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root:IsoMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_bcd_ratio'),
        (loc + '/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root:IsoMu24_PtEtaBins/efficienciesDATA/pt_abseta_DATA', 'm_trg_bcd_data_eff'),
        (loc + '/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root:IsoMu24_PtEtaBins/efficienciesMC/pt_abseta_MC', 'm_trg_mc_eff'),
        (loc + '/RunABCD_SF_ID.root:NUM_MediumID_DEN_TrackerMuons_pt_abseta', 'm_id_ratio'),
        (loc + '/RunABCD_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_pt_abseta', 'm_iso_ratio'),
    ]
    histsToWrapNewStyle = []

for task in histsToWrap:
    SafeWrapHist(w, ['m_pt', 'expr::m_abs_eta("TMath::Abs(@0)", m_eta[0])'], GetFromTFile(task[0]), name=task[1])

for task in histsToWrapNewStyle:
    SafeWrapHist(w, ['m_eta', 'm_pt'], GetFromTFile(task[0]), name=task[1])

if era == '2016':
    # Factor of 0.439 is the lumi of G+H
    for t in ['trg', 'id', 'iso']:
        w.factory('expr::m_%s_ratio("@0*(1-0.439) + @1*0.439", m_%s_bf_ratio, m_%s_gh_ratio)' % (t, t, t))

    for t in ['trg']:
        w.factory('expr::m_%s_data_eff("@0*(1-0.439) + @1*0.439", m_%s_bf_data_eff, m_%s_gh_data_eff)' % (t, t, t))

if era == '2018':
    #Lumi up to HLT update was 8.95, total 59.74
    for t in ['trg']:
        w.factory('expr::m_%s_ratio("@0*(0.1498)+@1*(1-0.1498)", m_%s_a_ratio, m_%s_bcd_ratio)' % (t, t, t))
        w.factory('expr::m_%s_data_eff("@0*(0.1498)+@1*(1-0.1498)", m_%s_a_data_eff, m_%s_bcd_data_eff)' % (t, t, t))

w.factory('expr::m_idiso_ratio("@0*@1", m_id_ratio, m_iso_ratio)')


###############################################################################
## Electrons
###############################################################################
loc = 'hvm/inputs/electrons/%s' % era

if era == '2016':
    histsToWrap = [
        (loc + '/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root:EGamma_SF2D', 'e_gsf_ratio'),
        (loc + '/2016LegacyReReco_ElectronMVA80_Fall17V2.root:EGamma_SF2D', 'e_id_ratio')
    ]
    histsToWrapInv = [
        (loc + '/ZEETP_2016_Data_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_data_eff'),
        (loc + '/ZEETP_2016_DY_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_mc_eff')
    ]
if era == '2017':
    histsToWrap = [
        (loc + '/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root:EGamma_SF2D', 'e_gsf_ratio'),
        (loc + '/2017_ElectronMVA80.root:EGamma_SF2D', 'e_id_ratio')
    ]
    histsToWrapInv = [
        (loc + '/ZEETP_2017_Data_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_data_eff'),
        (loc + '/ZEETP_2017_DY_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_mc_eff')
    ]

if era == '2018':
    histsToWrap = [
        (loc + '/egammaEffi.txt_EGM2D_updatedAll.root:EGamma_SF2D', 'e_gsf_ratio'),
        (loc + '/2018_ElectronMVA80.root:EGamma_SF2D', 'e_id_ratio')
    ]
    histsToWrapInv = [
        (loc + '/ZEETP_2018_Data_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_data_eff'),
        (loc + '/ZEETP_2018_DY_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_mc_eff')
    ]


for task in histsToWrap:
    SafeWrapHist(w, ['e_eta', 'e_pt'], GetFromTFile(task[0]), name=task[1])

for task in histsToWrapInv:
    SafeWrapHist(w, ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])'], GetFromTFile(task[0]), name=task[1])

w.factory('expr::e_gsfidiso_ratio("@0*@1", e_gsf_ratio, e_id_ratio)')
w.factory('expr::e_trg_ratio("@0/@1", e_trg_data_eff, e_trg_mc_eff)')

loc='hvm/inputs/rhoiso/%s' %era

if era == '2016':
    histsToWrap = [
       (loc + '/ZMMTP_2016_Data_Fits_TkIso_pt_bins_inc_eta.root:TkIso_pt_bins_inc_eta', 'rho_iso_data_eff_etainc'), 
       (loc + '/ZMMTP_2016_DY_Fits_TkIso_pt_bins_inc_eta.root:TkIso_pt_bins_inc_eta', 'rho_iso_mc_eff_etainc'),
       (loc + '/ZMMTP_2016_Data_Fits_TkIso_pt_eta_bins.root:TkIso_pt_eta_bins', 'rho_iso_data_eff'),
       (loc + '/ZMMTP_2016_DY_Fits_TkIso_pt_eta_bins.root:TkIso_pt_eta_bins', 'rho_iso_mc_eff')
    ]
if era == '2017':
    histsToWrap = [
       (loc + '/ZMMTP_2017_Data_Fits_TkIso_pt_bins_inc_eta.root:TkIso_pt_bins_inc_eta', 'rho_iso_data_eff_etainc'),
       (loc + '/ZMMTP_2017_DY_Fits_TkIso_pt_bins_inc_eta.root:TkIso_pt_bins_inc_eta', 'rho_iso_mc_eff_etainc'),
       (loc + '/ZMMTP_2017_Data_Fits_TkIso_pt_eta_bins.root:TkIso_pt_eta_bins', 'rho_iso_data_eff'),
       (loc + '/ZMMTP_2017_DY_Fits_TkIso_pt_eta_bins.root:TkIso_pt_eta_bins', 'rho_iso_mc_eff')
    ]
if era == '2018':
    histsToWrap = [
       (loc + '/ZMMTP_2018_Data_Fits_TkIso_pt_bins_inc_eta.root:TkIso_pt_bins_inc_eta', 'rho_iso_data_eff_etainc'),
       (loc + '/ZMMTP_2018_DY_Fits_TkIso_pt_bins_inc_eta.root:TkIso_pt_bins_inc_eta', 'rho_iso_mc_eff_etainc'),
       (loc + '/ZMMTP_2018_Data_Fits_TkIso_pt_eta_bins.root:TkIso_pt_eta_bins', 'rho_iso_data_eff'),
       (loc + '/ZMMTP_2018_DY_Fits_TkIso_pt_eta_bins.root:TkIso_pt_eta_bins', 'rho_iso_mc_eff')
    ]

for task in histsToWrap:
   SafeWrapHist(w, ['rho_pt', 'expr::rho_abs_eta("TMath::Abs(@0)",rho_eta[0])'], GetFromTFile(task[0]), name=task[1])

w.factory('expr::rhoiso_ratio_etainc("@0/@1", rho_iso_data_eff_etainc, rho_iso_mc_eff_etainc)')
w.factory('expr::rhoiso_ratio("@0/@1", rho_iso_data_eff, rho_iso_mc_eff)')

w.Print()
w.writeToFile('hvm_corrections_%s_v5.root' % era)
w.Delete()

