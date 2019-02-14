import ROOT
import json
# import sys
from pprint import pprint
from collections import defaultdict
import argparse
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)


parser = argparse.ArgumentParser()

parser.add_argument('--task', default='eft_region', choices=['eft_region', 'baseline', 'photon_fakes'])
parser.add_argument('--year', default='2016', choices=['2016', '2017', '2018'])
parser.add_argument('--indir', default='output/130818/wgamma_2016_v2/WGamma/')

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

year = args.year
doEle = True

remaps = {
    "2016": {
        'DY': 'DYJetsToLL_M-50-madgraphMLM',
        'ZG': 'ZGToLLG-amcatnloFXFX',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'SingleElectron',
        'TT': 'TT-powheg',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'W': 'WJetsToLNu-madgraphMLM',
        'VVTo2L2Nu': 'VVTo2L2Nu-amcatnloFXFX',
        'WWTo1L1Nu2Q': 'WWTo1L1Nu2Q-amcatnloFXFX',
        'WZTo1L1Nu2Q': 'WZTo1L1Nu2Q-amcatnloFXFX',
        'WZTo1L3Nu': 'WZTo1L3Nu-amcatnloFXFX',
        'WZTo2L2Q': 'WZTo2L2Q-amcatnloFXFX',
        'WZTo3LNu': 'WZTo3LNu-amcatnloFXFX',
        'ZZTo2L2Q': 'ZZTo2L2Q-amcatnloFXFX',
        'ZZTo4L': 'ZZTo4L-powheg',
        'ST_s': 'ST_s-channel_4f-amcatnlo',
        'ST_t_antitop': 'ST_t-channel_antitop_4f-powhegV2',
        'ST_t_top': 'ST_t-channel_top_4f-powhegV2',
        'ST_tW_antitop': 'ST_tW_antitop_5f-powheg',
        'ST_tW_top': 'ST_tW_top_5f-powheg',
        'TTG_DL': 'TTGamma_Dilept-madgraph',
        'TTG_Had': 'TTGamma_Hadronic-madgraph',
        'TTG_SL_T': 'TTGamma_SingleLeptFromT-madgraph',
        'TTG_SL_Tbar': 'TTGamma_SingleLeptFromTbar-madgraph',
    },
    "2017": {
        'DY': 'DYJetsToLL_M-50-madgraphMLM',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'SingleElectron',
        'TT_SL': 'TTToSemiLeptonic-powheg',
        'TT_Had': 'TTToHadronic-powheg',
        'TT_DL': 'TTTo2L2Nu-powheg',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'W': 'WJetsToLNu-madgraphMLM',
        'VVTo2L2Nu': 'VVTo2L2Nu-amcatnloFXFX',
        'WWTo1L1Nu2Q': 'WWTo1L1Nu2Q-amcatnloFXFX',
        'WZTo1L1Nu2Q': 'WZTo1L1Nu2Q-amcatnloFXFX',
        'WZTo1L3Nu': 'WZTo1L3Nu-amcatnloFXFX',
        'WZTo2L2Q': 'WZTo2L2Q-amcatnloFXFX',
        'WZTo3LNu': 'WZTo3LNu-amcatnloFXFX',
        'ZZTo2L2Q': 'ZZTo2L2Q-amcatnloFXFX',
        'ZZTo4L': 'ZZTo4L-powheg',
        'ST_s': 'ST_s-channel_4f-amcatnlo',
        'ST_t_antitop': 'ST_t-channel_antitop_4f-powhegV2',
        'ST_t_top': 'ST_t-channel_top_4f-powhegV2',
        'ST_tW_antitop': 'ST_tW_antitop_5f-powheg',
        'ST_tW_top': 'ST_tW_top_5f-powheg',
        'TTG_DL': 'TTGamma_Dilept-madgraph',
        'TTG_Had': 'TTGamma_Hadronic-madgraph',
        'TTG_SL_T': 'TTGamma_SingleLeptFromT-madgraph',
        'TTG_SL_Tbar': 'TTGamma_SingleLeptFromTbar-madgraph',
    },
    "2018": {
        'DY': 'DYJetsToLL_M-50-madgraphMLM',
        'ZG': 'ZGToLLG-amcatnloFXFX',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'EGamma',
        'TT_SL': 'TTToSemiLeptonic-powheg',
        'TT_Had': 'TTToHadronic-powheg',
        'TT_DL': 'TTTo2L2Nu-powheg',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'W': 'WJetsToLNu-madgraphMLM',
        'WW': 'WW-pythia',
        'WZ': 'WZ-pythia',
        'ZZ': 'ZZ-pythia',
        'ST_s': 'ST_s-channel_4f-amcatnlo',
        'ST_t_antitop': 'ST_t-channel_antitop_4f-powheg',
        'ST_t_top': 'ST_t-channel_top_4f-powheg',
        'ST_tW_antitop': 'ST_tW_antitop_5f-powheg',
        'ST_tW_top': 'ST_tW_top_5f-powheg',
        'TTG_DL': 'TTGamma_Dilept-madgraph',
        'TTG_Had': 'TTGamma_Hadronic-madgraph',
        'TTG_SL_T': 'TTGamma_SingleLeptFromT-madgraph',
        'TTG_SL_Tbar': 'TTGamma_SingleLeptFromTbar-madgraph',
    }
}

remap = remaps[year]

ttgammas = ['TTG_DL', 'TTG_Had', 'TTG_SL_T', 'TTG_SL_Tbar']
if year == '2016':
    dibosons = ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
    tt_samples = ['TT']
    zg_samples = ['ZG']
if year == '2017':
    dibosons = ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L']
    tt_samples = ['TT_SL', 'TT_Had', 'TT_DL']
    zg_samples = []
if year == '2018':
    dibosons = ['WW', 'WZ', 'ZZ']
    tt_samples = ['TT_SL', 'TT_Had', 'TT_DL']
    zg_samples = ['ZG']

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')

hists = Node()
do_cats = {}
do_cats['e'] = []
do_cats['m'] = []

X = SelectionManager()

if args.task == 'eft_region':
    pt_bins_min = [150, 210, 300, 420, 600, 850]
    pt_bins_max = [210, 300, 420, 600, 850, 1200]
    phi_bins_min = [0.00, 0.63, 1.26, 1.89, 2.52]
    phi_bins_max = [0.63, 1.26, 1.89, 2.52, 3.15]

    X['iso_t'] = 'p0_chiso < 0.441'
    X['iso_l'] = 'p0_chiso > 2 && p0_chiso < 10'
    X['sig_t'] = '(abs(p0_eta) < 1.4442 && p0_sigma < 0.01022) || (abs(p0_eta) > 1.4442 && p0_sigma < 0.03001)'
    X['sig_l'] = '(abs(p0_eta) < 1.4442 && p0_sigma > 0.01022) || (abs(p0_eta) > 1.4442 && p0_sigma > 0.03001)'
    X['baseline'] ='l0_pdgid == 13 && l0_trg && n_pre_m==1 && l0_iso<0.15 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && !p0_haspix && l0_pt>80 && met>80 && p0_pt>150 && l0p0_dr>3.0'
    if doEle:
        X['baseline'] ='l0_pdgid == 11 && l0_trg && n_pre_e==1 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && !p0_haspix && l0_pt>80 && met>80 && p0_pt>150 && l0p0_dr>3.0'
    X['baseline_wt'] = 'wt_def*wt_pu*wt_m0*wt_trg_m0*wt_p0'
    do_cats.extend(['baseline'])

    X['p_gen_ooa'] = 'gen_l0_q==+1 && !(gen_l0_pt>80 && gen_met>40 && gen_p0_pt>150 && gen_l0p0_dr>3.0)'
    X['n_gen_ooa'] = 'gen_l0_q==-1 && !(gen_l0_pt>80 && gen_met>40 && gen_p0_pt>150 && gen_l0p0_dr>3.0)'

    X['p_gen_acc'] = 'gen_l0_q==+1 && gen_l0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'
    X['n_gen_acc'] = 'gen_l0_q==-1 && gen_l0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'

    X['p_gen_acc_met1'] = 'gen_l0_q==+1 && gen_l0_pt>80 && gen_met>40 && gen_met<80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'
    X['n_gen_acc_met1'] = 'gen_l0_q==-1 && gen_l0_pt>80 && gen_met>40 && gen_met<80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'

    for i in range(len(pt_bins_min)):
        # Reconstructed categories
        X['p_%i' % i] = '$baseline && l0_q==+1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_%i' % i] = '$baseline && l0_q==-1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        do_cats.extend(['p_%i' % i, 'n_%i' % i])

        # Gen level selections
        X['p_gen_%i' % (i)] = '$p_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_gen_%i' % (i)] = '$n_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['p_gen_met1_%i' % (i)] = '$p_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_gen_met1_%i' % (i)] = '$n_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])

        for j in range(len(phi_bins_min)):
            X['p_gen_%i_%i' % (i, j)] = '$p_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f && abs(true_phi) >= %f && abs(true_phi) < %f' % (pt_bins_min[i], pt_bins_max[i], phi_bins_min[j], phi_bins_max[j])
            X['n_gen_%i_%i' % (i, j)] = '$n_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f && abs(true_phi) >= %f && abs(true_phi) < %f' % (pt_bins_min[i], pt_bins_max[i], phi_bins_min[j], phi_bins_max[j])
            X['p_gen_met1_%i_%i' % (i, j)] = '$p_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f && abs(true_phi) >= %f && abs(true_phi) < %f' % (pt_bins_min[i], pt_bins_max[i], phi_bins_min[j], phi_bins_max[j])
            X['n_gen_met1_%i_%i' % (i, j)] = '$n_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f && abs(true_phi) >= %f && abs(true_phi) < %f' % (pt_bins_min[i], pt_bins_max[i], phi_bins_min[j], phi_bins_max[j])

    drawvars = [
        ('p0_pt', pt_bins_min + [pt_bins_max[-1]]),
        ('abs(reco_phi)', (5, 0, 3.15)),
        ('abs(gen_phi)', (5, 0, 3.15)),
        ('abs(true_phi)', (5, 0, 3.15)),
    ]

    for sel in do_cats:
        for var, binning in drawvars:
            for sample in samples:
                hists[sel][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
            for P in ['W', 'DY'] + zg_samples + tt_samples + dibosons + ttgammas:
                hists[sel][var]['%s_R' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && p0_isprompt'), wt=X.get('$baseline_wt'))
                hists[sel][var]['%s_F' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && !p0_isprompt'), wt=X.get('$baseline_wt'))
            hists[sel][var]['data_fakes'] = Hist('TH1D', sample='data_obs', var=[var], binning=binning,
                sel=X.get('$' + sel, override={"sig_t": "$sig_l"}, printlevel=2),
                wt=X.get('$baseline_wt * wt_p0_fake'))
            hists[sel][var]['WG_p_ooa'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_ooa'), wt=X.get('$baseline_wt'))
            hists[sel][var]['WG_n_ooa'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_ooa'), wt=X.get('$baseline_wt'))
            hists[sel][var]['WG_p_acc'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_acc'), wt=X.get('$baseline_wt'))
            hists[sel][var]['WG_n_acc'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_acc'), wt=X.get('$baseline_wt'))
            hists[sel][var]['WG_p_acc_met1'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_acc_met1'), wt=X.get('$baseline_wt'))
            hists[sel][var]['WG_n_acc_met1'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_acc_met1'), wt=X.get('$baseline_wt'))
            if sel[-1].isdigit():
                # for i in range(len(pt_bins_min)):
                for i in [int(sel[-1])]:
                    hists[sel][var]['WG_p_%i' % (i)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_%i' % (i)), wt=X.get('$baseline_wt'))
                    hists[sel][var]['WG_n_%i' % (i)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_%i' % (i)), wt=X.get('$baseline_wt'))
                    hists[sel][var]['WG_p_met1_%i' % (i)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_met1_%i' % (i)), wt=X.get('$baseline_wt'))
                    hists[sel][var]['WG_n_met1_%i' % (i)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_met1_%i' % (i)), wt=X.get('$baseline_wt'))
                    for j in range(len(phi_bins_min)):
                        hists[sel][var]['WG_p_%i_%i' % (i, j)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_%i_%i' % (i, j)), wt=X.get('$baseline_wt'))
                        hists[sel][var]['WG_n_%i_%i' % (i, j)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_%i_%i' % (i, j)), wt=X.get('$baseline_wt'))
                        hists[sel][var]['WG_p_met1_%i_%i' % (i, j)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_met1_%i_%i' % (i, j)), wt=X.get('$baseline_wt'))
                        hists[sel][var]['WG_n_met1_%i_%i' % (i, j)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_met1_%i_%i' % (i, j)), wt=X.get('$baseline_wt'))

if args.task == 'baseline':
    X['iso_t'] = 'p0_chiso < 0.441'
    X['iso_l'] = 'p0_chiso > 2 && p0_chiso < 10'
    X['sig_t'] = '(abs(p0_eta) < 1.4442 && p0_sigma < 0.01022) || (abs(p0_eta) > 1.4442 && p0_sigma < 0.03001)'
    X['sig_l'] = '(abs(p0_eta) < 1.4442 && p0_sigma > 0.01022) || (abs(p0_eta) > 1.4442 && p0_sigma > 0.03001)'

    X['baseline_m_nopix'] ='l0_pdgid == 13 && l0_trg && n_pre_m==1 && l0_iso<0.15 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>30 && met>0 && p0_pt>30 && l0p0_dr>0.7 && l0met_mt>60'
    X['baseline_e_nopix'] ='l0_pdgid == 11 && l0_trg && n_pre_e==1 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>35 && met>0 && p0_pt>30 && l0p0_dr>0.7 && l0met_mt>60'
    X['baseline_m'] ='$baseline_m_nopix && !p0_haspix'
    X['baseline_e'] ='$baseline_e_nopix && !p0_haspix'
    X['baseline_wt'] = 'wt_def*wt_pu*wt_m0*wt_trg_m0*wt_p0'
    X['baseline_M_veto_m'] ='$baseline_m && abs(91.2 - l0p0_M) > 15'
    X['baseline_M_veto_e'] ='$baseline_e && abs(91.2 - l0p0_M) > 15'

    do_cats['e'].extend(['baseline_e_nopix', 'baseline_e', 'baseline_M_veto_e'])
    do_cats['m'].extend(['baseline_m_nopix', 'baseline_m', 'baseline_M_veto_m'])

    drawvars = [
        ('n_vtx', (30, 0., 60.)),
        ('l0met_mt', (30, 0., 200.)),
        ('l0_pt', (40, 0., 150.)),
        ('l0_eta', (20, -3.0, 3.0)),
        ('l0_phi', (20, -3.15, 3.15)),
        ('l0_iso', (40, 0, 2.0)),
        ('l1_pt', (40, 0., 150.)),
        ('l1_eta', (20, -3.0, 3.0)),
        ('l0l1_M', (40, 60, 120)),
        ('l0l1_dr', (20, 0., 5.)),
        ('met', (20, 0., 200.)),
        ('met_phi', (20, -3.15, 3.15)),
        ('xy_met', (20, 0., 200.)),
        ('xy_met_phi', (20, -3.15, 3.15)),
        ('puppi_met', (20, 0., 200.)),
        ('puppi_met_phi', (20, -3.15, 3.15)),
        ('p0_pt', [0, 10, 20, 30, 40, 50, 60, 80, 100, 120, 160, 200, 250, 300]),
        ('p0_eta', (20, -3.0, 3.0)),
        ('p0_phi', (30, -3.15, 3.15)),
        ('l0p0_dr', (20, 0., 5.)),
        ('l0p0_M', (20, 60, 120)),
        ('p0_chiso', (40, 0, 20.0)),
        ('p0_neiso', (40, 0, 20.0)),
        ('p0_phiso', (40, 0, 20.0)),
        ('p0_hovere', (20, 0., 0.5)),
        ('p0_sigma', (30, 0., 0.06)),
        ('p0_haspix', (2, -0.5, 1.5)),
        ('wt_def', (100, 0, 2)),
        ('wt_pu', (100, 0, 2)),
        ('wt_m0', (100, 0, 2)),
        ('wt_trg_m0', (100, 0, 2)),
        ('wt_p0', (100, 0, 2))
    ]

    for chn in ['e', 'm']:
        for sel in do_cats[chn]:
            for var, binning in drawvars:
                for sample in samples:
                    hists[chn][sel][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                for P in ['W', 'DY'] + zg_samples + tt_samples + dibosons + ttgammas:
                    hists[chn][sel][var]['%s_R' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && p0_isprompt'), wt=X.get('$baseline_wt'))
                    hists[chn][sel][var]['%s_F' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && !p0_isprompt'), wt=X.get('$baseline_wt'))
                hists[chn][sel][var]['data_obs'] = Hist('TH1D', sample='data_obs_%s' % chn, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                hists[chn][sel][var]['data_fakes'] = Hist('TH1D', sample='data_obs_%s' % chn, var=[var], binning=binning,
                    sel=X.get('$' + sel, override={"sig_t": "$sig_l"}),
                    wt=X.get('$baseline_wt * wt_p0_fake'))
            for sample in samples:
                hists['2D'][chn][sel]['l0_eta_phi'][sample] = Hist('TH2F', sample=sample, var=['l0_eta', 'l0_phi'], binning=(30, -3, 3, 30, -3.15, 3.15), sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                hists['2D'][chn][sel]['p0_eta_phi'][sample] = Hist('TH2F', sample=sample, var=['p0_eta', 'p0_phi'], binning=(30, -3, 3, 30, -3.15, 3.15), sel=X.get('$' + sel), wt=X.get('$baseline_wt'))


if args.task == 'photon_fakes':
    X['baseline'] = 'm0_trg && n_m==1 && m0_iso<0.15 && m0met_mt>60 && n_p==1 && m0p0_dr>0.7 && p0_medium_noch && !p0_haspix'
    X['baseline_wt'] = 'wt_def*wt_pu*wt_m0*wt_trg_m0*wt_p0'
    X['barrel'] = '$baseline && abs(p0_eta) < 1.4442'
    X['endcap'] = '$baseline && abs(p0_eta) > 1.4442'
    X['iso_t'] = 'p0_chiso < 0.441'
    X['iso_l'] = 'p0_chiso > 2 && p0_chiso < 10'
    X['sig_t'] = '(abs(p0_eta) < 1.4442 && p0_sigma < 0.01022) || (abs(p0_eta) > 1.4442 && p0_sigma < 0.03001)'
    X['sig_l'] = '(abs(p0_eta) < 1.4442 && p0_sigma > 0.01022) || (abs(p0_eta) > 1.4442 && p0_sigma > 0.03001)'
    do_cats.extend(['baseline', 'barrel', 'endcap'])
    for S in ['barrel', 'endcap']:
        X['%s_iso_l' % S] = '$%s && $iso_l' % S
        X['%s_sig_l' % S] = '$%s && $sig_l' % S
        X['%s_iso_l_sig_t' % S] = '$%s && $iso_l && $sig_t' % S
        X['%s_iso_l_sig_l' % S] = '$%s && $iso_l && $sig_l' % S
        X['%s_iso_t_sig_l' % S] = '$%s && $iso_t && $sig_l' % S
        X['%s_iso_t_sig_t' % S] = '$%s && $iso_t && $sig_t' % S
        do_cats.extend(['%s_iso_l' % S, '%s_sig_l' % S, '%s_iso_l_sig_t' % S, '%s_iso_l_sig_l' % S, '%s_iso_t_sig_l' % S, '%s_iso_t_sig_t' % S])

    drawvars = [
        ('m0met_mt', (30, 0., 200.)),
        ('m0_pt', (40, 0., 150.)),
        ('m0_eta', (20, -3.0, 3.0)),
        ('m0_iso', (40, 0, 2.0)),
        ('m1_pt', (40, 0., 150.)),
        ('m1_eta', (20, -3.0, 3.0)),
        ('m0m1_M', (40, 60, 120)),
        ('m0m1_dr', (20, 0., 5.)),
        ('met', (20, 0., 200.)),
        ('p0_pt', (100, 0, 500.)),
        ('p0_eta', (20, -3.0, 3.0)),
        ('m0p0_dr', (20, 0., 5.)),
        ('m0p0_M', (20, 60, 120)),
        ('p0_chiso', (40, 0, 20.0)),
        ('p0_neiso', (40, 0, 20.0)),
        ('p0_phiso', (40, 0, 20.0)),
        ('p0_hovere', (20, 0., 0.5)),
        ('p0_sigma', (30, 0., 0.06)),
        ('p0_haspix', (2, -0.5, 1.5))
    ]

    for sel in do_cats:
        for var, binning in drawvars:
            for sample in samples:
                hists[sel][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                for P in ['W', 'DY'] + zg_samples + tt_samples + dibosons + ttgammas:
                    hists[sel][var]['%s_R' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && p0_isprompt'), wt=X.get('$baseline_wt'))
                    hists[sel][var]['%s_F' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && !p0_isprompt'), wt=X.get('$baseline_wt'))
                if 'iso_t_sig_t' in sel:
                    hists[sel][var]['data_fakes'] = Hist('TH1D', sample='data_obs', var=[var], binning=binning,
                        sel=X.get('$' + sel, override={"sig_t": "$sig_l"}, printlevel=2),
                        wt=X.get('$baseline_wt * wt_p0_fake'))


MultiDraw(hists, samples, tname, mt_cores=4)

with open('input/cfg_wgamma_%s_v3.json' % year) as jsonfile:
    cfg = json.load(jsonfile)
    sample_cfg = cfg['samples']

for sample in samples:
    f = ROOT.TFile.Open(samples[sample])
    sample_cfg[remap[sample]]['events'] = f.Get('counters').GetBinContent(2)
    f.Close()

for chn in ['e', 'm']:
    for path, hname, obj in hists[chn].ListObjects():
        name = obj.sample
        if 'data_obs' not in name:
            tgt_lumi = sample_cfg[remap['data_obs_%s' % chn]]['lumi']
            events = sample_cfg[remap[name]]['events']
            xsec = sample_cfg[remap[name]]['xsec']
            scale = tgt_lumi * xsec / events
            obj.Scale(scale)
        if 'data_obs' in name:
            tgt_lumi = sample_cfg[remap['data_obs_%s' % chn]]['lumi']
            obj.SetTitle('lumi:%.1f fb^{-1} (13 TeV, %s)' % (tgt_lumi / 1000., year))
    for sample in samples:
        if 'data_obs' not in sample:
            tgt_lumi = sample_cfg[remap['data_obs_%s' % chn]]['lumi']
            events = sample_cfg[remap[sample]]['events']
            xsec = sample_cfg[remap[sample]]['xsec']
            scale = tgt_lumi * xsec / events
            print '%-50s %-20.1f %-12.2f %-12.2f' % (remap[sample], events, xsec, scale)


fout = ROOT.TFile('output_%s_%s.root' % (year, args.task), 'RECREATE')

for path, node in hists.ListNodes(withObjects=True):
    print path
    if path.startswith('2D'):
        continue
    node['VV'] = sum([node[X] for X in dibosons[1:]], node[dibosons[0]])
    node['VV_F'] = sum([node['%s_F' % X] for X in dibosons[1:]], node['%s_F' % dibosons[0]])
    node['VV_R'] = sum([node['%s_R' % X] for X in dibosons[1:]], node['%s_R' % dibosons[0]])

    node['TTG'] = sum([node[X] for X in ttgammas[1:]], node[ttgammas[0]])
    node['TTG_F'] = sum([node['%s_F' % X] for X in ttgammas[1:]], node['%s_F' % ttgammas[0]])
    node['TTG_R'] = sum([node['%s_R' % X] for X in ttgammas[1:]], node['%s_R' % ttgammas[0]])

    if year in ['2017', '2018']:
        node['TT'] = sum([node[X] for X in tt_samples[1:]], node[tt_samples[0]])
        node['TT_F'] = sum([node['%s_F' % X] for X in tt_samples[1:]], node['%s_F' % tt_samples[0]])
        node['TT_R'] = sum([node['%s_R' % X] for X in tt_samples[1:]], node['%s_R' % tt_samples[0]])

    node['Total_R'] = node['WG'] + node['TT_R'] + node['DY_R'] + node['VV_R']
    node['Total_F'] = node['W_F'] + node['TT_F'] + node['DY_F'] + node['VV_F']

NodeToTDir(fout, hists)

fout.Close()
