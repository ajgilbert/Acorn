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

parser.add_argument('--task', default='eft_region', choices=['eft_region', 'baseline', 'photon_fakes', 'electron_fakes'])
parser.add_argument('--label', default='default')
parser.add_argument('--year', default='2016', choices=['2016', '2017', '2018'])
parser.add_argument('--indir', default='output/130818/wgamma_2016_v2/WGamma/')
parser.add_argument('--extra-cfg', default=None)

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

year = args.year

remaps = {
    "2016": {
        'DY': 'DYJetsToLL_M-50-madgraphMLM',
        'ZG': 'ZGToLLG-amcatnloFXFX',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'SingleElectron',
        'TT': 'TT-powheg',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG_NLO': 'WGToLNuG-amcatnloFXFX-stitched',
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
        'GG': 'DiPhotonJetsBox_MGG-80toInf'
    },
    "2017": {
        'DY': 'DYJetsToLL_M-50-madgraphMLM',
        'ZG': 'ZGToLLG-amcatnloFXFX',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'SingleElectron',
        'TT_SL': 'TTToSemiLeptonic-powheg',
        'TT_Had': 'TTToHadronic-powheg',
        'TT_DL': 'TTTo2L2Nu-powheg',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG_NLO': 'WGToLNuG-amcatnloFXFX-stitched',
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
        'GG': 'DiPhotonJetsBox_MGG-80toInf'
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
        'WG_NLO': 'WGToLNuG-amcatnloFXFX-stitched',
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
        'GG': 'DiPhotonJetsBox_MGG-80toInf'
    }
}

remap = remaps[year]

wg_sample = 'WG_NLO'
wg_samples = ['WG', 'WG_NLO']
w_samples = ['W']
ttg_samples = ['TTG_DL', 'TTG_Had', 'TTG_SL_T', 'TTG_SL_Tbar']
zg_samples = ['ZG']
dy_samples = ['DY']

if year == '2016':
    vv_samples = ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
    tt_samples = ['TT']
if year == '2017':
    vv_samples = ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
    tt_samples = ['TT_SL', 'TT_Had', 'TT_DL']
if year == '2018':
    vv_samples = ['WW', 'WZ', 'ZZ', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
    tt_samples = ['TT_SL', 'TT_Had', 'TT_DL']

all_samples = wg_samples + w_samples + ttg_samples + tt_samples + zg_samples + dy_samples + vv_samples

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')


def AddPhotonSplitting(node, label, sample, var_list, binning, sel, wt, components=['R', 'J', 'E']):
    if 'R' in components:
        node['%s_R' % label] = Hist('TH1F', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && $p0_isprompt' % sel), wt=X.get(wt))
    if 'J' in components:
        node['%s_J' % label] = Hist('TH1F', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && $p0_isjet' % sel), wt=X.get(wt))
    if 'E' in components:
        node['%s_E' % label] = Hist('TH1F', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && $p0_iselec' % sel), wt=X.get(wt))


def AddDYSplitting(node, label, sample, var_list, binning, sel, wt, components=['R', 'J', 'E'], postfix=''):
    node[label + '_XZG' + postfix] = Hist('TH1F', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && !(gen_is_zg && p0_truth==1)' % sel), wt=X.get(wt))
    node[label + '_IZG' + postfix] = Hist('TH1F', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && (gen_is_zg && p0_truth==1)' % sel), wt=X.get(wt))
    AddPhotonSplitting(node, label + '_XZG' + postfix, sample, var_list, binning, '(%s) && !(gen_is_zg && p0_truth==1)' % sel, wt, components=components)
    AddPhotonSplitting(node, label + '_IZG' + postfix, sample, var_list, binning, '(%s) && (gen_is_zg && p0_truth==1)' % sel, wt, components=components)


def HistSum(label_list):
    return sum([node[X] for X in label_list[1:]], node[label_list[0]])


def CapNegativeBins(h):
    for ibin in xrange(1, h.GetNbinsX() + 1):
        if h.GetBinContent(ibin) < 0.:
            h.SetBinContent(ibin, 0.)


hists = Node()
do_cats = {}
do_cats['e'] = []
do_cats['m'] = []

X = SelectionManager()
# Selections for photon truth
X['p0_isprompt'] = 'p0_truth == 1 || p0_truth == 4 || p0_truth == 5'
X['p0_isjet'] = 'p0_truth == 6 || p0_truth == 0 || p0_truth == 3'
X['p0_iselec'] = 'p0_truth == 2'
X['p0_isfake'] = 'p0_truth == 6 || p0_truth == 0 || p0_truth == 2 || p0_truth == 3'

# Selecting barrel or endcap photons
X['p0_eb'] = 'abs(p0_eta) < 1.4442'
X['p0_ee'] = 'abs(p0_eta) > 1.4442'

# Selections for photon isolation and sigma_ietaieta
X['iso_t'] = '($p0_eb && p0_chiso < 1.141) || ($p0_ee && p0_chiso < 1.051)'
X['iso_l'] = 'p0_chiso > 4 && p0_chiso < 13'
X['sig_t'] = '($p0_eb && p0_sigma < 0.01015) || ($p0_ee && p0_sigma < 0.0272)'
X['sig_l'] = '($p0_eb && p0_sigma > 0.01015) || ($p0_ee && p0_sigma > 0.0272)'

# Selections for vetoing e->photon fakes
X['efake_veto_e'] = '!p0_haspix && p0_eveto && n_veto_e == 1'
X['efake_veto_m'] = '!p0_haspix && p0_eveto && n_veto_m == 1'

# Selection to veto (or select) mZ region
X['mZ_veto'] = 'abs(91.2 - l0p0_M) > 15'
X['mZ_veto_inv'] = 'abs(91.2 - l0p0_M) <= 15'

# Analysis selection levels:
X['baseline_m_nomt_nopix'] ='l0_pdgid == 13 && l0_trg && n_pre_m==1 && l0_iso<0.15 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>30 && met>0 && p0_pt>30 && l0p0_dr>0.7'
X['baseline_e_nomt_nopix'] ='l0_pdgid == 11 && l0_trg && n_pre_e==1 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>35 && met>0 && p0_pt>30 && l0p0_dr>0.7'
X['baseline_m_nopix'] ='$baseline_m_nomt_nopix && l0met_mt>60'
X['baseline_e_nopix'] ='$baseline_e_nomt_nopix && l0met_mt>60'
X['baseline_m_nomt'] ='$baseline_m_nomt_nopix && $efake_veto_m'
X['baseline_e_nomt'] ='$baseline_e_nomt_nopix && $efake_veto_e'
X['baseline_m'] ='$baseline_m_nopix && $efake_veto_m'
X['baseline_e'] ='$baseline_e_nopix && $efake_veto_e'
X['baseline_m_mZ_veto'] ='$baseline_m'
X['baseline_e_mZ_veto'] ='$baseline_e && $mZ_veto'
X['baseline_m_mZ_veto_t'] ='$baseline_m_mZ_veto && l0_tight'
X['baseline_e_mZ_veto_t'] ='$baseline_e_mZ_veto && l0_tight'

# EFT region
X['baseline_m_eft'] ='l0_pdgid == 13 && l0_trg && n_pre_m==1 && l0_iso<0.15 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>80 && met>80 && p0_pt>150 && l0p0_dr>3.0 && $efake_veto_m'
X['baseline_e_eft'] ='l0_pdgid == 11 && l0_trg && n_pre_e==1 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>80 && met>80 && p0_pt>150 && l0p0_dr>3.0 && $mZ_veto && $efake_veto_e'

# Event weights
X['baseline_wt'] = 'wt_def*wt_pu*wt_l0*wt_trg_l0*wt_p0*wt_p0_e_fake*wt_pf'


if args.task == 'eft_region':
    eft_defaults = {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        # 'phi_var': 'abs(true_phi)',
        # 'phi_bins': '(5,0.,math.pi)'
        'phi_var': 'abs(gen_sphi)',
        'phi_var_obs': 'abs(reco_sphi)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(6,0.,math.pi/2.)',
    }

    if args.extra_cfg is not None:
        eft_defaults.update(json.loads(args.extra_cfg))
    pprint(eft_defaults)
    pt_bins = BinEdgesFromStr(eft_defaults['pt_bins'])
    phi_bins = BinEdgesFromStr(eft_defaults['phi_bins'])
    pt_bins_min = pt_bins[:-1]
    pt_bins_max = pt_bins[1:]
    phi_bins_min = phi_bins[:-1]
    phi_bins_max = phi_bins[1:]
    phi_var = eft_defaults['phi_var']

    do_cats['e'].extend(['baseline_e_eft'])
    do_cats['m'].extend(['baseline_m_eft'])

    X['p_gen_ooa'] = 'gen_l0_q==+1 && !(gen_l0_pt>80 && gen_met>40 && gen_p0_pt>150 && gen_l0p0_dr>3.0)'
    X['n_gen_ooa'] = 'gen_l0_q==-1 && !(gen_l0_pt>80 && gen_met>40 && gen_p0_pt>150 && gen_l0p0_dr>3.0)'

    X['p_gen_acc'] = 'gen_l0_q==+1 && gen_l0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'
    X['n_gen_acc'] = 'gen_l0_q==-1 && gen_l0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'

    X['p_gen_acc_met1'] = 'gen_l0_q==+1 && gen_l0_pt>80 && gen_met>40 && gen_met<80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'
    X['n_gen_acc_met1'] = 'gen_l0_q==-1 && gen_l0_pt>80 && gen_met>40 && gen_met<80 && gen_p0_pt>150 && gen_l0p0_dr>3.0'

    for i in range(len(pt_bins_min)):
        # Reconstructed categories
        X['p_e_%i' % i] = '$baseline_e_eft && l0_q==+1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_e_%i' % i] = '$baseline_e_eft && l0_q==-1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['p_m_%i' % i] = '$baseline_m_eft && l0_q==+1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_m_%i' % i] = '$baseline_m_eft && l0_q==-1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        do_cats['e'].extend(['p_e_%i' % i, 'n_e_%i' % i])
        do_cats['m'].extend(['p_m_%i' % i, 'n_m_%i' % i])

        # Gen level selections
        X['p_gen_%i' % (i)] = '$p_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_gen_%i' % (i)] = '$n_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['p_gen_met1_%i' % (i)] = '$p_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_gen_met1_%i' % (i)] = '$n_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])

        for j in range(len(phi_bins_min)):
            X['p_gen_%i_%i' % (i, j)] = '$p_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['n_gen_%i_%i' % (i, j)] = '$n_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['p_gen_met1_%i_%i' % (i, j)] = '$p_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['n_gen_met1_%i_%i' % (i, j)] = '$n_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])

    drawvars = [
        ('gen_p0_pt', BinningFromStr(eft_defaults['pt_bins'])),
        ('p0_pt', BinningFromStr(eft_defaults['pt_bins'])),
        (eft_defaults['phi_var'], BinningFromStr(eft_defaults['phi_bins'])),
        (eft_defaults['phi_var_obs'], BinningFromStr(eft_defaults['phi_bins_obs']))
    ]

    for chn in ['e', 'm']:
        for sel in do_cats[chn]:
            for var, binning in drawvars:
                for sample in samples:
                    hists[chn][sel][var][sample] = Hist('TH1F', sample=sample, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                for P in all_samples:
                    AddPhotonSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')
                    AddPhotonSplitting(hists[chn][sel][var], P + '_fw', P, [var], binning, X.get('$' + sel, override={"sig_t": "$sig_l"}), '$baseline_wt * wt_p0_fake', components=['R', 'E'])
                for P in dy_samples + zg_samples:
                    AddDYSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')
                    AddDYSplitting(hists[chn][sel][var], P, P, [var], binning, X.get('$' + sel, override={"sig_t": "$sig_l"}), '$baseline_wt * wt_p0_fake', components=['R', 'E'], postfix='_fw')
                hists[chn][sel][var]['data_obs'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                hists[chn][sel][var]['data_fakes'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning,
                    sel=X.get('$' + sel, override={"sig_t": "$sig_l"}, printlevel=0),
                    wt=X.get('$baseline_wt * wt_p0_fake'))
                for P in [wg_sample]:
                    hists[chn][sel][var]['WG_ooa_p'] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_ooa'), wt=X.get('$baseline_wt'))
                    hists[chn][sel][var]['WG_ooa_n'] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_ooa'), wt=X.get('$baseline_wt'))
                    hists[chn][sel][var]['WG_main_p'] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_acc'), wt=X.get('$baseline_wt'))
                    hists[chn][sel][var]['WG_main_n'] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_acc'), wt=X.get('$baseline_wt'))
                    hists[chn][sel][var]['WG_met1_p'] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $p_gen_acc_met1'), wt=X.get('$baseline_wt'))
                    hists[chn][sel][var]['WG_met1_n'] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $n_gen_acc_met1'), wt=X.get('$baseline_wt'))
                    if sel[-1].isdigit():
                        chg = sel[0]
                        for i in range(len(pt_bins_min)):
                            hists[chn][sel][var]['WG_main_%s_%i' % (chg, i)] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_%i' % (chg, i)), wt=X.get('$baseline_wt'))
                            hists[chn][sel][var]['WG_met1_%s_%i' % (chg, i)] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_met1_%i' % (chg, i)), wt=X.get('$baseline_wt'))
                        for i in [int(sel[-1])]:
                            for j in range(len(phi_bins_min)):
                                hists[chn][sel][var]['WG_main_%s_%i_%i' % (chg, i, j)] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_%i_%i' % (chg, i, j)), wt=X.get('$baseline_wt'))
                                hists[chn][sel][var]['WG_met1_%s_%i_%i' % (chg, i, j)] = Hist('TH1F', sample=P, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_met1_%i_%i' % (chg, i, j)), wt=X.get('$baseline_wt'))
        pdgid = '11' if chn == 'e' else '13'
        for var, binning in drawvars:
            hists[chn]['XS'][var]['XS_WG_p_%s_acc' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $p_gen_acc' % pdgid), wt=X.get('wt_def'))
            hists[chn]['XS'][var]['XS_WG_n_%s_acc' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $n_gen_acc' % pdgid), wt=X.get('wt_def'))
        hists[chn]['XS']['2D']['XS_WG_p_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $p_gen_acc' % pdgid), wt=X.get('wt_def'))
        hists[chn]['XS']['2D']['XS_WG_n_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $n_gen_acc' % pdgid), wt=X.get('wt_def'))
            # # Pick this up automatically
            # hists[chn]['XS'][var]['XS_WG_p_%s_acc_sphi_0' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $p_gen_acc && abs(gen_sphi)>=0 && abs(gen_sphi)<(TMath::Pi()/6)' % pdgid), wt=X.get('wt_def'))
            # hists[chn]['XS'][var]['XS_WG_n_%s_acc_sphi_0' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $n_gen_acc && abs(gen_sphi)>=0 && abs(gen_sphi)<(TMath::Pi()/6)' % pdgid), wt=X.get('wt_def'))
            # hists[chn]['XS'][var]['XS_WG_p_%s_acc_sphi_1' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $p_gen_acc && abs(gen_sphi)>=(TMath::Pi()/6) && abs(gen_sphi)<(TMath::Pi()/3)' % pdgid), wt=X.get('wt_def'))
            # hists[chn]['XS'][var]['XS_WG_n_%s_acc_sphi_1' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $n_gen_acc && abs(gen_sphi)>=(TMath::Pi()/6) && abs(gen_sphi)<(TMath::Pi()/3)' % pdgid), wt=X.get('wt_def'))
            # hists[chn]['XS'][var]['XS_WG_p_%s_acc_sphi_2' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $p_gen_acc && abs(gen_sphi)>=(TMath::Pi()/3) && abs(gen_sphi)<(TMath::Pi()/2)' % pdgid), wt=X.get('wt_def'))
            # hists[chn]['XS'][var]['XS_WG_n_%s_acc_sphi_2' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $n_gen_acc && abs(gen_sphi)>=(TMath::Pi()/3) && abs(gen_sphi)<(TMath::Pi()/2)' % pdgid), wt=X.get('wt_def'))

if args.task in ['baseline', 'electron_fakes']:
    if args.task == 'electron_fakes':
        pt_bins_min = [30, 35, 40, 60, 100]
        pt_bins_max = [35, 40, 60, 100, 200]
        eta_bins_min = [0, 1.4442]
        eta_bins_max = [1.4442, 2.5]
        X['baseline_wt'] = 'wt_def*wt_pu*wt_l0*wt_trg_l0*wt_p0*wt_pf'

        for e_min, e_max in zip(eta_bins_min, eta_bins_max):
            for p_min, p_max in zip(pt_bins_min, pt_bins_max):
                label = '%g_%g_%g_%g' % (e_min, e_max, p_min, p_max)
                label = label.replace('.', 'p')
                X['fail_%s' % label] ='$baseline_e_nopix && $mZ_veto_inv && !$efake_veto_e && abs(p0_eta) >= %g && abs(p0_eta) < %g && p0_pt >= %g && p0_pt < %g' % (e_min, e_max, p_min, p_max)
                X['pass_%s' % label] ='$baseline_e_nopix && $mZ_veto_inv && $efake_veto_e && abs(p0_eta) >= %g && abs(p0_eta) < %g && p0_pt >= %g && p0_pt < %g' % (e_min, e_max, p_min, p_max)
                do_cats['e'].append('fail_%s' % label)
                do_cats['e'].append('pass_%s' % label)

    do_cats['e'].extend(['baseline_e_nomt_nopix', 'baseline_e_nomt', 'baseline_e_nopix', 'baseline_e', 'baseline_e_mZ_veto'])
    do_cats['m'].extend(['baseline_m_nomt_nopix', 'baseline_m_nomt', 'baseline_m_nopix', 'baseline_m', 'baseline_m_mZ_veto'])

    drawvars = [
        ('n_vtx', (30, 0., 60.)),
        ('n_veto_m', (3, -0.5, 2.5)),
        ('n_veto_e', (3, -0.5, 2.5)),
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
        ('p0_eveto', (2, -0.5, 1.5)),
        ('p0_truth', (7, -0.5, 6.5)),
        ('wt_def', (100, 0, 2)),
        ('wt_pu', (100, 0, 2)),
        ('wt_l0', (100, 0, 2)),
        ('wt_trg_l0', (100, 0, 2)),
        ('wt_p0', (100, 0, 2)),
        ('wt_p0_e_fake', (100, 0, 4))
    ]
    if args.task == 'electron_fakes':
        drawvars = [
            ('l0p0_M', (40, 60, 120)),
        ]

    for chn in ['e', 'm']:
        for sel in do_cats[chn]:
            for var, binning in drawvars:
                for sample in samples:
                    hists[chn][sel][var][sample] = Hist('TH1F', sample=sample, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                for P in all_samples:
                    AddPhotonSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')
                    AddPhotonSplitting(hists[chn][sel][var], P + '_fw', P, [var], binning, X.get('$' + sel, override={"sig_t": "$sig_l"}), '$baseline_wt * wt_p0_fake', components=['R', 'E'])
                for P in dy_samples + zg_samples:
                    AddDYSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')
                    AddDYSplitting(hists[chn][sel][var], P, P, [var], binning, X.get('$' + sel, override={"sig_t": "$sig_l"}), '$baseline_wt * wt_p0_fake', components=['R', 'E'], postfix='_fw')

                hists[chn][sel][var]['data_obs'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                hists[chn][sel][var]['data_fakes'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning,
                    sel=X.get('$' + sel, override={"sig_t": "$sig_l"}),
                    wt=X.get('$baseline_wt * wt_p0_fake'))
            for sample in samples:
                hists['2D'][chn][sel]['l0_eta_phi'][sample] = Hist('TH2F', sample=sample, var=['l0_eta', 'l0_phi'], binning=(30, -3, 3, 30, -3.15, 3.15), sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                hists['2D'][chn][sel]['p0_eta_phi'][sample] = Hist('TH2F', sample=sample, var=['p0_eta', 'p0_phi'], binning=(30, -3, 3, 30, -3.15, 3.15), sel=X.get('$' + sel), wt=X.get('$baseline_wt'))


if args.task == 'photon_fakes':
    X['baseline_m_mZ_veto'] = X.get('$baseline_m_mZ_veto', override={"sig_t": "1", "iso_t": "1"})
    X['baseline_e_mZ_veto'] = X.get('$baseline_e_mZ_veto', override={"sig_t": "1", "iso_t": "1"})
    X['barrel_m'] = '$baseline_m_mZ_veto && $p0_eb'
    X['endcap_m'] = '$baseline_m_mZ_veto && $p0_ee'
    X['barrel_e'] = '$baseline_e_mZ_veto && $p0_eb'
    X['endcap_e'] = '$baseline_e_mZ_veto && $p0_ee'

    # do_cats['e'].extend(['baseline_e_mZ_veto', 'barrel_e', 'endcap_e'])
    # do_cats['m'].extend(['baseline_m_mZ_veto', 'barrel_m', 'endcap_m'])

    for chn in ['e', 'm']:
        for S in ['barrel_%s' % chn, 'endcap_%s' % chn]:
            X['%s_iso_l' % S] = '$%s && $iso_l' % S
            X['%s_sig_l' % S] = '$%s && $sig_l' % S
            X['%s_iso_l_sig_t' % S] = '$%s && $iso_l && $sig_t' % S
            X['%s_iso_l_sig_l' % S] = '$%s && $iso_l && $sig_l' % S
            X['%s_iso_t_sig_l' % S] = '$%s && $iso_t && $sig_l' % S
            X['%s_iso_t_sig_t' % S] = '$%s && $iso_t && $sig_t' % S
            # X['%s_p0_medium' % S] = '$%s && p0_medium' % S
            do_cats[chn].extend(['%s_iso_l' % S, '%s_sig_l' % S, '%s_iso_l_sig_t' % S, '%s_iso_l_sig_l' % S, '%s_iso_t_sig_l' % S, '%s_iso_t_sig_t' % S])

    drawvars = [
        ('l0met_mt', (30, 0., 200.)),
        ('l0_pt', (40, 0., 150.)),
        ('l0_eta', (20, -3.0, 3.0)),
        # ('l0l1_M', (40, 60, 120)),
        ('met', (20, 0., 200.)),
        ('p0_pt', (100, 0, 500.)),
        ('p0_eta', (20, -3.0, 3.0)),
        # ('l0p0_dr', (20, 0., 5.)),
        # ('l0p0_M', (20, 60, 120)),
        ('p0_chiso', (40, 0, 20.0)),
        # ('p0_neiso', (40, 0, 20.0)),
        # ('p0_phiso', (40, 0, 20.0)),
        # ('p0_hovere', (20, 0., 0.5)),
        ('p0_sigma', (30, 0., 0.06)),
        ('p0_haspix', (2, -0.5, 1.5)),
        ('p0_truth', (7, -0.5, 6.5)),
    ]

    for chn in ['e', 'm']:
        for sel in do_cats[chn]:
            for var, binning in drawvars:
                for sample in samples:
                    hists[chn][sel][var][sample] = Hist('TH1F', sample=sample, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                for P in all_samples:
                    AddPhotonSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')
                    if 'iso_t_sig_t' in sel:
                        AddPhotonSplitting(hists[chn][sel][var], P + '_fw', P, [var], binning, X.get('$' + sel, override={"sig_t": "$sig_l"}), '$baseline_wt * wt_p0_fake', components=['R', 'E'])
                for P in dy_samples + zg_samples:
                    AddDYSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')
                    if 'iso_t_sig_t' in sel:
                        AddDYSplitting(hists[chn][sel][var], P, P, [var], binning, X.get('$' + sel, override={"sig_t": "$sig_l"}), '$baseline_wt * wt_p0_fake', components=['R', 'E'], postfix='_fw')
                    hists[chn][sel][var]['data_obs'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                    if 'iso_t_sig_t' in sel:
                        hists[chn][sel][var]['data_fakes'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning,
                            sel=X.get('$' + sel, override={"sig_t": "$sig_l"}, printlevel=0),
                            wt=X.get('$baseline_wt * wt_p0_fake'))


MultiDraw(hists, samples, tname, mt_cores=4, mt_thresh=1)

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
            if hname.startswith('XS_'):
                print path, hname
                tgt_lumi = 1.0
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


fout = ROOT.TFile('output_%s_%s_%s.root' % (year, args.task, args.label), 'RECREATE')

for path, node in hists.ListNodes(withObjects=True):
    print path
    if path.startswith('2D'):
        continue
    if 'XS' in path.split('/'):
        continue

    node['VV'] = HistSum(vv_samples)
    node['TTG'] = HistSum(ttg_samples)
    if year in ['2017', '2018']:
        node['TT'] = HistSum(tt_samples)

    for P in ['E', 'J', 'R']:
        node['VV_%s' % P] = HistSum(['%s_%s' % (lbl, P) for lbl in vv_samples])
        node['TTG_%s' % P] = HistSum(['%s_%s' % (lbl, P) for lbl in ttg_samples])
        if year in ['2017', '2018']:
            node['TT_%s' % P] = HistSum(['%s_%s' % (lbl, P) for lbl in tt_samples])

    node['Total_R'] = HistSum([wg_sample, 'TT_R', 'DY_XZG_R', 'ZG_IZG_R', 'VV_R'])
    node['Total_E'] = HistSum(['W_E', 'TT_E', 'DY_E', 'VV_E'])
    node['Total_J'] = HistSum(['W_J', 'TT_J', 'DY_J', 'VV_J'])

    if 'data_fakes' in node.d:
        all_r_samples = [wg_sample, 'DY_XZG', 'ZG_IZG'] + tt_samples + vv_samples
        all_e_samples = dy_samples + w_samples + tt_samples + vv_samples
        node['Total_fw_R'] = HistSum([lbl + '_fw_R' for lbl in all_r_samples])
        node['Total_fw_E'] = HistSum([lbl + '_fw_E' for lbl in all_e_samples])
        node['data_fakes_sub'] = node['data_fakes'] - (node['Total_fw_R'] + node['Total_fw_E'])
        CapNegativeBins(node['data_fakes_sub'])

    for S in ['Total', 'W', 'TT', 'DY', 'VV']:
        node['%s_F' % S] = HistSum(['%s_E' % S, '%s_J' % S])

NodeToTDir(fout, hists)

fout.Close()
