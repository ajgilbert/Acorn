#!/usr/bin/env python
import ROOT
# import json
import argparse
from array import array
from Acorn.Analysis.workspaceTools import *
import CombineHarvester.CombineTools.plotting as plot

plot.ModTDRStyle()

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
    mc_edges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
                34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
                66, 67, 68, 69, 70, 71, 72, 73, 74]
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
    h_data = GetFromTFile('wgamma/inputs/pileup_wgamma_2016_v3.root:pileup')
    h_data_hi = GetFromTFile('wgamma/inputs/pileup_wgamma_2016_v3_hi.root:pileup')
    h_data_lo = GetFromTFile('wgamma/inputs/pileup_wgamma_2016_v3_lo.root:pileup')

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
    h_data = GetFromTFile('wgamma/inputs/pileup_wgamma_2017_v3.root:pileup')
    h_data_hi = GetFromTFile('wgamma/inputs/pileup_wgamma_2017_v3_hi.root:pileup')
    h_data_lo = GetFromTFile('wgamma/inputs/pileup_wgamma_2017_v3_lo.root:pileup')

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
    h_data = GetFromTFile('wgamma/inputs/pileup_wgamma_2018_v4.root:pileup')
    h_data_hi = GetFromTFile('wgamma/inputs/pileup_wgamma_2018_v4_hi.root:pileup')
    h_data_lo = GetFromTFile('wgamma/inputs/pileup_wgamma_2018_v4_lo.root:pileup')

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
h_data_hi.Scale(1. / h_data_hi.Integral())
h_data_lo.Scale(1. / h_data_lo.Integral())
h_mc.Scale(1. / h_mc.Integral())

SafeWrapHist(w, ['pu_int'], h_data, name='pileup_data')
SafeWrapHist(w, ['pu_int'], h_data_hi, name='pileup_data_hi')
SafeWrapHist(w, ['pu_int'], h_data_lo, name='pileup_data_lo')
SafeWrapHist(w, ['pu_int'], h_mc, name='pileup_mc')

w.factory('expr::pileup_ratio("@0/@1", pileup_data, pileup_mc)')
w.factory('expr::pileup_ratio_hi("@0/@1", pileup_data_hi, pileup_mc)')
w.factory('expr::pileup_ratio_lo("@0/@1", pileup_data_lo, pileup_mc)')

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
loc = 'wgamma/inputs/muons/%s' % era
wrapvars_def = ['m_pt', 'expr::m_abs_eta("TMath::Abs(@0)",m_eta[0])']
wrapvars_new = ['m_eta', 'm_pt']
if era == '2016':
    histsToWrap = [
        (loc + '/trigger/EfficienciesAndSF_RunBtoF.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_bf_ratio', wrapvars_def),
        (loc + '/trigger/EfficienciesAndSF_Period4.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_gh_ratio', wrapvars_def),
        (loc + '/trigger/EfficienciesAndSF_RunBtoF.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_bf_ratio_err', wrapvars_def, +1.0),
        (loc + '/trigger/EfficienciesAndSF_Period4.root:IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_gh_ratio_err', wrapvars_def, +1.0),
        (loc + '/id/RunBCDEF_SF_ID.root:NUM_MediumID_DEN_genTracks_eta_pt', 'm_id_bf_ratio', wrapvars_new),
        (loc + '/id/RunGH_SF_ID.root:NUM_MediumID_DEN_genTracks_eta_pt', 'm_id_gh_ratio', wrapvars_new),
        (loc + '/id/RunBCDEF_SF_ID.root:NUM_MediumID_DEN_genTracks_eta_pt', 'm_id_bf_ratio_err', wrapvars_new, +1.0),
        (loc + '/id/RunGH_SF_ID.root:NUM_MediumID_DEN_genTracks_eta_pt', 'm_id_gh_ratio_err', wrapvars_new, +1.0),
        (loc + '/iso/RunBCDEF_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_eta_pt', 'm_iso_bf_ratio', wrapvars_new),
        (loc + '/iso/RunGH_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_eta_pt', 'm_iso_gh_ratio', wrapvars_new),
        (loc + '/iso/RunBCDEF_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_eta_pt', 'm_iso_bf_ratio_err', wrapvars_new, +1.0),
        (loc + '/iso/RunGH_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_eta_pt', 'm_iso_gh_ratio_err', wrapvars_new, +1.0),
        (loc + '/output_2016_lepton_fakes_ratios_m_200511.root:lepton_fakes', 'm_fake_ratio', wrapvars_def + ['l0met_mt']),
    ]

    muon_trk_eff_hist = TGraphAsymmErrorsToTH1D(GetFromTFile(loc + '/track/Tracking_EfficienciesAndSF_BCDEFGH.root:ratio_eff_eta3_dr030e030_corr'))
    SafeWrapHist(w, ['m_eta'], muon_trk_eff_hist, name='m_trk_ratio')

if era == '2017':
    histsToWrap = [
        (loc + '/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root:IsoMu27_PtEtaBins/pt_abseta_ratio', 'm_trg_ratio', wrapvars_def),
        (loc + '/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root:IsoMu27_PtEtaBins/pt_abseta_ratio', 'm_trg_ratio_err', wrapvars_def, +1.0),
        (loc + '/RunBCDEF_SF_ID_syst.root:NUM_MediumID_DEN_genTracks_pt_abseta', 'm_id_ratio', wrapvars_def),
        (loc + '/RunBCDEF_SF_ID_syst.root:NUM_MediumID_DEN_genTracks_pt_abseta', 'm_id_ratio_err', wrapvars_def, +1.0),
        (loc + '/RunBCDEF_SF_ISO_syst.root:NUM_TightRelIso_DEN_MediumID_pt_abseta', 'm_iso_ratio', wrapvars_def),
        (loc + '/RunBCDEF_SF_ISO_syst.root:NUM_TightRelIso_DEN_MediumID_pt_abseta', 'm_iso_ratio_err', wrapvars_def, +1.0),
        (loc + '/output_2017_lepton_fakes_ratios_m_200511.root:lepton_fakes', 'm_fake_ratio', wrapvars_def + ['l0met_mt']),
    ]

if era == '2018':
    histsToWrap = [
        (loc + '/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root:NUM_MediumID_DEN_TrackerMuons_pt_abseta', 'm_id_ratio', wrapvars_def),
        (loc + '/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root:NUM_MediumID_DEN_TrackerMuons_pt_abseta', 'm_id_ratio_err', wrapvars_def, +1.0),
        (loc + '/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_pt_abseta', 'm_iso_ratio', wrapvars_def),
        (loc + '/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ISO.root:NUM_TightRelIso_DEN_MediumID_pt_abseta', 'm_iso_ratio_err', wrapvars_def, +1.0),
        (loc + '/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root:IsoMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_early_ratio', wrapvars_def),
        (loc + '/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root:IsoMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_early_ratio_err', wrapvars_def, +1.0),
        (loc + '/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root:IsoMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_late_ratio', wrapvars_def),
        (loc + '/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root:IsoMu24_PtEtaBins/pt_abseta_ratio', 'm_trg_late_ratio_err', wrapvars_def, +1.0),
        (loc + '/output_2018_lepton_fakes_ratios_m_200511.root:lepton_fakes', 'm_fake_ratio', wrapvars_def + ['l0met_mt']),
    ]

for task in histsToWrap:
    if len(task) > 3:
        SafeWrapHist(w, task[2], HistErr(GetFromTFile(task[0]), task[3]), name=task[1])
    else:
        SafeWrapHist(w, task[2], GetFromTFile(task[0]), name=task[1])


if era == '2016':
    # Factor of 0.439 is the lumi of G+H
    for t in ['trg', 'id', 'iso']:
        w.factory('expr::m_%s_ratio("@0*(1-0.439) + @1*0.439", m_%s_bf_ratio, m_%s_gh_ratio)' % (t, t, t))
        # This formula is for 100% uncorrelated between run periods
        # w.factory('expr::m_%s_ratio_err("TMath::Sqrt(@0*@0*@1*@1*0.561*0.561 + @2*@2*@3*@3*0.439*0.439)/@4", m_%s_bf_ratio, m_%s_bf_ratio_err, m_%s_gh_ratio, m_%s_gh_ratio_err, m_%s_ratio)' % (t, t, t, t, t, t))
        # This just averages between the two errors
        w.factory('expr::m_%s_ratio_err("@0*(1-0.439) + @1*0.439", m_%s_bf_ratio_err, m_%s_gh_ratio_err)' % (t, t, t))

if era == '2017':
    w.factory('expr::m_trk_ratio("1", m_eta[0])')

if era == '2018':
    w.factory('expr::m_trk_ratio("1", m_eta[0])')
    w.factory('expr::m_trg_ratio("@0*(1-0.85) + @1*0.85", m_trg_early_ratio, m_trg_late_ratio)')
    w.factory('expr::m_trg_ratio_err("TMath::Sqrt(@0*@0*@1*@1*0.15*0.15 + @2*@2*@3*@3*0.85*0.85)/@4", m_trg_early_ratio, m_trg_early_ratio_err, m_trg_late_ratio, m_trg_late_ratio_err, m_trg_ratio)')


w.factory('expr::m_idisotrk_ratio("@0*@1*@2", m_id_ratio, m_iso_ratio, m_trk_ratio)')
w.factory('expr::m_idisotrk_ratio_err("TMath::Sqrt(@0*@0+@1*@1)", m_id_ratio_err, m_iso_ratio_err)')


###############################################################################
## Electrons
###############################################################################
loc = 'wgamma/inputs/electrons/%s' % era

if era == '2016':
    histsToWrap = [
        (loc + '/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root:EGamma_SF2D', 'e_gsf_ratio', ['e_eta', 'e_pt']),
        (loc + '/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root:EGamma_SF2D', 'e_gsf_ratio_err', ['e_eta', 'e_pt'], +1.0),
        (loc + '/2016LegacyReReco_ElectronMedium.root:EGamma_SF2D', 'e_id_ratio', ['e_eta', 'e_pt']),
        (loc + '/2016LegacyReReco_ElectronMedium.root:EGamma_SF2D', 'e_id_ratio_err', ['e_eta', 'e_pt'], +1.0),
        (loc + '/ZeeTP_2016_Data_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_data', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])']),
        (loc + '/ZeeTP_2016_DYJetsToLL_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_mc', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])']),
        (loc + '/output_2016_lepton_fakes_ratios_e_200511.root:lepton_fakes', 'e_fake_ratio', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])', 'l0met_mt']),
    ]
if era == '2017':
    histsToWrap = [
        (loc + '/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root:EGamma_SF2D', 'e_gsf_ratio', ['e_eta', 'e_pt']),
        (loc + '/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root:EGamma_SF2D', 'e_gsf_ratio_err', ['e_eta', 'e_pt'], +1.0),
        (loc + '/2017_ElectronMedium.root:EGamma_SF2D', 'e_id_ratio', ['e_eta', 'e_pt']),
        (loc + '/2017_ElectronMedium.root:EGamma_SF2D', 'e_id_ratio_err', ['e_eta', 'e_pt'], +1.0),
        (loc + '/ZeeTP_2017_Data_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_data', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])']),
        (loc + '/ZeeTP_2017_DYJetsToLL_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_mc', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])']),
        (loc + '/output_2017_lepton_fakes_ratios_e_200511.root:lepton_fakes', 'e_fake_ratio', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])', 'l0met_mt']),
    ]

if era == '2018':
    histsToWrap = [
        (loc + '/egammaEffi.txt_EGM2D.root:EGamma_SF2D', 'e_gsf_ratio', ['e_eta', 'e_pt']),
        (loc + '/egammaEffi.txt_EGM2D.root:EGamma_SF2D', 'e_gsf_ratio_err', ['e_eta', 'e_pt'], +1.0),
        (loc + '/2018_ElectronMedium.root:EGamma_SF2D', 'e_id_ratio', ['e_eta', 'e_pt']),
        (loc + '/2018_ElectronMedium.root:EGamma_SF2D', 'e_id_ratio_err', ['e_eta', 'e_pt'], +1.0),
        (loc + '/ZeeTP_2018_Data_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_data', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])']),
        (loc + '/ZeeTP_2018_DYJetsToLL_Fits_Trg_pt_eta_bins.root:Trg_pt_eta_bins', 'e_trg_mc', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])']),
        (loc + '/output_2018_lepton_fakes_ratios_e_200511.root:lepton_fakes', 'e_fake_ratio', ['e_pt', 'expr::e_abs_eta("TMath::Abs(@0)",e_eta[0])', 'l0met_mt']),
    ]

for task in histsToWrap:
    if len(task) > 3:
        SafeWrapHist(w, task[2], HistErr(GetFromTFile(task[0]), task[3]), name=task[1])
    else:
        SafeWrapHist(w, task[2], GetFromTFile(task[0]), name=task[1])

w.factory('expr::e_gsfidiso_ratio("@0*@1", e_gsf_ratio, e_id_ratio)')
w.factory('expr::e_gsfidiso_ratio_err("TMath::Sqrt(@0*@0+@1*@1)", e_gsf_ratio_err, e_id_ratio_err)')
w.factory('expr::e_trg_ratio("@0/@1", e_trg_data, e_trg_mc)')

###############################################################################
## Photons
###############################################################################
loc = 'wgamma/inputs/photons/%s' % era

if era == '2016':
    histsToWrap = [
        (loc + '/Fall17V2_2016_Medium_photons.root:EGamma_SF2D', 'p_id_ratio', ['p_eta', 'p_pt']),
        (loc + '/Fall17V2_2016_Medium_photons.root:EGamma_EffData2D', 'p_id_data', ['p_eta', 'p_pt']),
        (loc + '/Fall17V2_2016_Medium_photons.root:EGamma_EffMC2D', 'p_id_mc', ['p_eta', 'p_pt']),
        (loc + '/Fall17V2_2016_Medium_photons.root:EGamma_SF2D', 'p_id_ratio_err', ['p_eta', 'p_pt'], +1.0),
        (loc + '/PhotonTP_2016_Data_Fits_WorstIso_pt_eta_bins.root:WorstIso_pt_eta_bins', 'p_swi_data', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        (loc + '/PhotonTP_2016_DYJetsToLL_Fits_WorstIso_pt_eta_bins.root:WorstIso_pt_eta_bins', 'p_swi_mc', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        # (loc + '/EFakesTP_2016_data_obs_Fits_EGammaFakes.root:EGammaFakes', 'e_p_fake_data', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        # (loc + '/EFakesTP_2016_DY_E_Fits_EGammaFakes.root:EGammaFakes', 'e_p_fake_mc', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])'])
    ]
    w.factory('expr::p_psv_ratio("(TMath::Abs(@0) < 1.4442) * 0.9978 + (TMath::Abs(@0) >= 1.4442) * 0.9931", p_eta[0], p_pt[0])')

if era == '2017':
    histsToWrap = [
        (loc + '/2017_PhotonsMedium.root:EGamma_SF2D', 'p_id_ratio', ['p_eta', 'p_pt']),
        (loc + '/2017_PhotonsMedium.root:EGamma_EffData2D', 'p_id_data', ['p_eta', 'p_pt']),
        (loc + '/2017_PhotonsMedium.root:EGamma_EffMC2D', 'p_id_mc', ['p_eta', 'p_pt']),
        (loc + '/2017_PhotonsMedium.root:EGamma_SF2D', 'p_id_ratio_err', ['p_eta', 'p_pt'], +1.0),
        (loc + '/PhotonTP_2017_Data_Fits_WorstIso_pt_eta_bins.root:WorstIso_pt_eta_bins', 'p_swi_data', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        (loc + '/PhotonTP_2017_DYJetsToLL_Fits_WorstIso_pt_eta_bins.root:WorstIso_pt_eta_bins', 'p_swi_mc', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        # (loc + '/EFakesTP_2017_data_obs_Fits_EGammaFakes.root:EGammaFakes', 'e_p_fake_data', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        # (loc + '/EFakesTP_2017_DY_E_Fits_EGammaFakes.root:EGammaFakes', 'e_p_fake_mc', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])'])
    ]
    w.factory('expr::p_psv_ratio("(TMath::Abs(@0) < 1.4442) * 0.96726 + (TMath::Abs(@0) >= 1.4442) * 0.914782", p_eta[0], p_pt[0])')

if era == '2018':
    histsToWrap = [
        (loc + '/2018_PhotonsMedium.root:EGamma_SF2D', 'p_id_ratio', ['p_eta', 'p_pt']),
        (loc + '/2018_PhotonsMedium.root:EGamma_EffData2D', 'p_id_data', ['p_eta', 'p_pt']),
        (loc + '/2018_PhotonsMedium.root:EGamma_EffMC2D', 'p_id_mc', ['p_eta', 'p_pt']),
        (loc + '/2018_PhotonsMedium.root:EGamma_SF2D', 'p_id_ratio_err', ['p_eta', 'p_pt'], +1.0),
        (loc + '/HasPix_2018.root:eleVeto_SF', 'p_psv_ratio', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        (loc + '/PhotonTP_2018_Data_Fits_WorstIso_pt_eta_bins.root:WorstIso_pt_eta_bins', 'p_swi_data', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        (loc + '/PhotonTP_2018_DYJetsToLL_Fits_WorstIso_pt_eta_bins.root:WorstIso_pt_eta_bins', 'p_swi_mc', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        # (loc + '/EFakesTP_2018_data_obs_nominal_Fits_EGammaFakes.root:EGammaFakes', 'e_p_fake_data', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])']),
        # (loc + '/EFakesTP_2018_DY_E_Fits_EGammaFakes.root:EGammaFakes', 'e_p_fake_mc', ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])'])

    ]

h_efake_nominal = GetFromTFile(loc + '/EFakesTP_%s_data_obs_nominal_Fits_EGammaFakes.root:EGammaFakes' % era)
h_efake_hi = GetFromTFile(loc + '/EFakesTP_%s_data_obs_bkgHi_Fits_EGammaFakes.root:EGammaFakes' % era)
h_efake_lo = GetFromTFile(loc + '/EFakesTP_%s_data_obs_bkgLo_Fits_EGammaFakes.root:EGammaFakes' % era)
h_efake_mc = GetFromTFile(loc + '/EFakesTP_%s_DY_E_nominal_Fits_EGammaFakes.root:EGammaFakes' % era)

SafeWrapHist(w, ['p_pt', 'p_eta'], h_efake_nominal, name='e_p_fake_data_stat')
SafeWrapHist(w, ['p_pt', 'p_eta'], h_efake_mc, name='e_p_fake_mc_stat')

h_efake_ratio_stat = h_efake_nominal.Clone()
h_efake_ratio_stat.Divide(h_efake_mc)
SafeWrapHist(w, ['p_pt', 'p_eta'], h_efake_ratio_stat, name='e_p_fake_ratio_stat')

h_efake_data_syst = HistSystVariations(h_efake_nominal, h_efake_lo, h_efake_hi, keepNominalErr=False)
SafeWrapHist(w, ['p_pt', 'p_eta'], h_efake_data_syst, name='e_p_fake_data_syst')

h_efake_ratio_syst = h_efake_data_syst.Clone()
h_efake_ratio_syst.Divide(ZeroErrors(h_efake_mc))
SafeWrapHist(w, ['p_pt', 'p_eta'], h_efake_ratio_syst, name='e_p_fake_ratio_syst')

h_efake_nominal = HistSystVariations(h_efake_nominal, h_efake_lo, h_efake_hi)
SafeWrapHist(w, ['p_pt', 'p_eta'], h_efake_nominal, name='e_p_fake_data')

h_efake_ratio = h_efake_nominal.Clone()
h_efake_ratio.Divide(h_efake_mc)
SafeWrapHist(w, ['p_pt', 'p_eta'], h_efake_ratio, name='e_p_fake_ratio')
SafeWrapHist(w, ['p_pt', 'p_eta'], HistErr(h_efake_ratio, +1.0), name='e_p_fake_ratio_err')

for task in histsToWrap:
    if len(task) > 3:
        SafeWrapHist(w, task[2], HistErr(GetFromTFile(task[0]), task[3]), name=task[1])
    else:
        SafeWrapHist(w, task[2], GetFromTFile(task[0]), name=task[1])


w.factory('expr::p_swi_ratio("@0/@1", p_swi_data, p_swi_mc)')

SummaryPlots({
    "h_ref": w.genobj('hist_p_id_ratio'),
    "proj": 'Y',
    "logx": False,
    "main_label": "p_id_%s" % era,
    "main_text": "Photon ID efficiency",
    "x_axis_title": "Photon candidate p_{T} (GeV)",
    "data": [w.genobj('hist_p_id_data')],
    "mc": [w.genobj('hist_p_id_mc')],
    "ratio": [w.genobj('hist_p_id_ratio')],
    "ratio_range": [0.7, 1.3],
    "y_range": [0, 0.2],
    "y_label": '#eta'
    })


SummaryPlots({
    "h_ref": h_efake_nominal,
    "proj": 'X',
    "logx": True,
    "main_label": "e_p_fake_%s" % era,
    "main_text": "Misid. e #rightarrow #gamma",
    "x_axis_title": "Photon candidate p_{T} (GeV)",
    "data": [w.genobj('hist_e_p_fake_data'), w.genobj('hist_e_p_fake_data_stat'), w.genobj('hist_e_p_fake_data_syst')],
    "mc": [w.genobj('hist_e_p_fake_mc_stat')],
    "ratio": [w.genobj('hist_e_p_fake_ratio'), w.genobj('hist_e_p_fake_ratio_stat'), w.genobj('hist_e_p_fake_ratio_syst')],
    "ratio_range": [0, 3.9],
    "y_range": [0, 0.2],
    "y_label": '#eta'
    })


SummaryPlots({
    "h_ref": h_efake_nominal,
    "proj": 'Y',
    "logx": False,
    "main_label": "e_p_fake_%s" % era,
    "main_text": "Misid. e #rightarrow #gamma",
    "x_axis_title": "Photon candidate #eta",
    "data": [w.genobj('hist_e_p_fake_data'), w.genobj('hist_e_p_fake_data_stat'), w.genobj('hist_e_p_fake_data_syst')],
    "mc": [w.genobj('hist_e_p_fake_mc_stat')],
    "ratio": [w.genobj('hist_e_p_fake_ratio'), w.genobj('hist_e_p_fake_ratio_stat'), w.genobj('hist_e_p_fake_ratio_syst')],
    "ratio_range": [0, 3.9],
    "y_range": [0, 0.2],
    "y_label": 'p_{T}'
    })

###############################################################################
## Photon fakes
###############################################################################
loc = 'wgamma/inputs/photons/%s' % era

# e and m fake weights have been merged now
histsToWrap = [
    (loc + '/output_%s_photon_fakes_ratios_m.root:photon_fakes_stat_syst' % era, 'p_fake_ratio'),
    (loc + '/output_%s_photon_fakes_ratios_m.root:photon_fakes_index' % era, 'p_fake_index'),
    (loc + '/output_%s_photon_fakes_ratios_m.root:photon_fakes_stat_syst' % era, 'p_fake_ratio_err', +1.0),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes' % era, 'p_fake_ratio_new'),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes_stat' % era, 'p_fake_ratio_stat_new'),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes_bkg_syst' % era, 'p_fake_ratio_bkg_syst_new'),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes_const_syst' % era, 'p_fake_ratio_const_syst_new'),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes_index' % era, 'p_fake_index_new'),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes' % era, 'p_fake_ratio_err_new', +1.0),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes_mc' % era, 'p_fake_ratio_mc_new'),
    (loc + '/output_%s_photon_fakes_ratios_SWI.root:photon_fakes_mc_true' % era, 'p_fake_ratio_mc_true_new'),
    (loc + '/high_pt_photon_fakes_MC.root:photon_fakes', 'p_highpt_fake_ratio'),
    (loc + '/high_pt_photon_fakes_MC.root:photon_fakes', 'p_highpt_fake_ratio_err', +1.0),
]

for task in histsToWrap:
    if len(task) > 2:
        SafeWrapHist(w, ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])'], HistErr(GetFromTFile(task[0]), task[2]), name=task[1])
    else:
        SafeWrapHist(w, ['p_pt', 'expr::p_abs_eta("TMath::Abs(@0)",p_eta[0])'], GetFromTFile(task[0]), name=task[1])

SummaryPlotsPhotonFakes({
    "h_ref": w.genobj('hist_p_fake_ratio_new'),
    "proj": 'X',
    "logx": False,
    "main_label": "p_fake_%s" % era,
    "main_text": "Misid. j #rightarrow #gamma",
    "x_axis_title": "Photon candidate p_{T} (GeV)",
    "data": [w.genobj('hist_p_fake_ratio_new'), w.genobj('hist_p_fake_ratio_stat_new'), w.genobj('hist_p_fake_ratio_const_syst_new'), w.genobj('hist_p_fake_ratio_bkg_syst_new')],
    "data_labels": ['total', 'stat', 'const_syst', 'bkg_syst'],
    "mc": [],
    "ratio": [],
    "ratio_range": [0, 3.9],
    "y_range": [0, 0.2],
    "y_label": 'p_{T}'
    })

SummaryPlotsPhotonFakes({
    "h_ref": w.genobj('hist_p_fake_ratio_mc_new'),
    "proj": 'X',
    "logx": False,
    "main_label": "p_fake_mc_%s" % era,
    "main_text": "Misid. j #rightarrow #gamma",
    "x_axis_title": "Photon candidate p_{T} (GeV)",
    "data": [w.genobj('hist_p_fake_ratio_mc_new'), w.genobj('hist_p_fake_ratio_mc_true_new')],
    "data_labels": ['mc', 'mc_true'],
    "mc": [],
    "ratio": [],
    "ratio_range": [0, 3.9],
    "y_range": [0, 0.2],
    "y_label": 'p_{T}'
    })

w.Print()
w.writeToFile('wgamma_corrections_%s_v12.root' % era)
w.Delete()

