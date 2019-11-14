import ROOT
import json
import sys
from pprint import pprint
from collections import defaultdict
import argparse
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
ROOT.TH1.SetDefaultSumw2()

parser = argparse.ArgumentParser()

parser.add_argument('--task', default='eft_region', choices=['eft_region', 'baseline', 'photon_fakes', 'photon_fakes2', 'electron_fakes', 'fid_region'])
parser.add_argument('--label', default='default')
parser.add_argument('--year', default='2016', choices=['2016', '2017', '2018'])
parser.add_argument('--indir', default='output/130818/wgamma_2016_v2/WGamma/')
parser.add_argument('--indir-data', default=None)
parser.add_argument('--syst', default=None)
parser.add_argument('--extra-cfg', default=None)

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

year = args.year

remaps = {
    "2016": {
        'DY': 'DYJetsToLL_M-50-amcatnloFXFX',
        'ZG': 'ZGToLLG-amcatnloFXFX',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'SingleElectron',
        'TT': 'TT-powheg',
        # 'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG': 'WGToLNuG-amcatnloFXFX-stitched',
        'W': 'WJetsToLNu-stitched',
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
        'GG': 'DiPhotonJetsBox_MGG-80toInf',
        'GG_NLO': 'DiPhotonJets_MGG-80toInf',
        'GG_L': 'DiPhotonJetsBox_MGG-40to80',
        'WWG': 'WWG-amcatnlo'
    },
    "2017": {
        'DY': 'DYJetsToLL_M-50-amcatnloFXFX',
        'ZG': 'ZGToLLG-amcatnloFXFX',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'SingleElectron',
        'TT_SL': 'TTToSemiLeptonic-powheg',
        'TT_Had': 'TTToHadronic-powheg',
        'TT_DL': 'TTTo2L2Nu-powheg',
        # 'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG': 'WGToLNuG-amcatnloFXFX-stitched',
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
        'GG_LO': 'DiPhotonJetsBox_MGG-80toInf',
        'GG_L': 'DiPhotonJetsBox_MGG-40to80'
    },
    "2018": {
        'DY': 'DYJetsToLL_M-50-amcatnloFXFX',
        'ZG': 'ZGToLLG-amcatnloFXFX',
        'data_obs_m': 'SingleMuon',
        'data_obs_e': 'EGamma',
        'TT_SL': 'TTToSemiLeptonic-powheg',
        'TT_Had': 'TTToHadronic-powheg',
        'TT_DL': 'TTTo2L2Nu-powheg',
        # 'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG': 'WGToLNuG-amcatnloFXFX-stitched',
        'W': 'WJetsToLNu-madgraphMLM',
        'VVTo2L2Nu': 'VVTo2L2Nu-amcatnloFXFX',
        'WWTo1L1Nu2Q': 'WWTo1L1Nu2Q-amcatnloFXFX',
        # 'WZTo1L1Nu2Q': 'WZTo1L1Nu2Q-amcatnloFXFX',
        'WZTo1L3Nu': 'WZTo1L3Nu-amcatnloFXFX',
        'WZTo2L2Q': 'WZTo2L2Q-amcatnloFXFX',
        'WZTo3LNu': 'WZTo3LNu-amcatnloFXFX',
        'ZZTo2L2Q': 'ZZTo2L2Q-amcatnloFXFX',
        'ZZTo4L': 'ZZTo4L-powheg',
        'ST_s': 'ST_s-channel_4f-amcatnlo',
        'ST_t_antitop': 'ST_t-channel_antitop_4f-powheg',
        'ST_t_top': 'ST_t-channel_top_4f-powheg',
        'ST_tW_antitop': 'ST_tW_antitop_5f-powheg',
        'ST_tW_top': 'ST_tW_top_5f-powheg',
        'TTG_DL': 'TTGamma_Dilept-madgraph',
        'TTG_Had': 'TTGamma_Hadronic-madgraph',
        'TTG_SL_T': 'TTGamma_SingleLeptFromT-madgraph',
        'TTG_SL_Tbar': 'TTGamma_SingleLeptFromTbar-madgraph',
        'GG': 'DiPhotonJetsBox_MGG-80toInf',
        'WWG': 'WWG-amcatnlo'
    }
}

remap = remaps[year]

wg_sample = 'WG'

# Combine samples together at the end
groups = {
    'WG': ['WG'],
    'W': ['W'],
    'DY': ['DY'],
    'ZG': ['ZG'],
    'TT': ['TT_SL', 'TT_Had', 'TT_DL'],
    'TTG': ['TTG_DL', 'TTG_Had', 'TTG_SL_T', 'TTG_SL_Tbar'],
    'VV': ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top'],
    'GG': ['GG_LO', 'GG_L'],
    'WWG': ['WWG']
}

if year == '2016':
    groups['TT'] = ['TT']
    groups['GG'] = ['GG_NLO', 'GG_L']

if year == '2018':
    groups['VV'] = ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
    groups['GG'] = ['GG']

samples = {}
for sa in remap:
    if 'data_obs' in sa and args.indir_data is not None:
        samples[sa] = (args.indir_data + remap[sa] + '.root')
    else:
        samples[sa] = (prefix + remap[sa] + '.root')

main_wt_systs = {
    'e': [
        ('wt_pf', 'CMS_prefiring'),
        ('wt_l0', 'CMS_eff_e'),
        ('wt_trg_l0', 'CMS_trigger_e'),
        ('wt_p0', 'CMS_eff_p'),
        ('wt_p0_e_fake', 'CMS_ele_fake_p'),
    ],
    'm': [
        ('wt_pf', 'CMS_prefiring'),
        ('wt_l0', 'CMS_eff_m'),
        ('wt_trg_l0', 'CMS_trigger_m'),
        ('wt_p0', 'CMS_eff_p'),
        ('wt_p0_e_fake', 'CMS_ele_fake_p'),
    ]
}


if args.syst is not None:
    # No weight-based systematics in this case
    main_wt_systs = {'e': [], 'm': []}

# Naming scheme for processes
# [MAIN]_[IZG]_[R/F/E]_[fw]_[CMS_eff_lUp]

def ApplyDYSplitting(hlist, components=['XZG', 'IZG']):
    res = []
    for label, sel, wt in hlist:
        if 'XZG' in components:
            res.append(('%s_XZG' % label, '(%s) && !(gen_is_zg && p0_truth==1)' % sel, wt))
        if 'IZG' in components:
            res.append(('%s_IZG' % label, '(%s) && (gen_is_zg && p0_truth==1)' % sel, wt))
    return res


def ApplyTTSplitting(hlist, components=['XTTG', 'ITTG']):
    res = []
    for label, sel, wt in hlist:
        if 'XTTG' in components:
            res.append(('%s_XTTG' % label, '(%s) && !(p0_truth==1)' % sel, wt))
        if 'ITTG' in components:
            res.append(('%s_ITTG' % label, '(%s) && (p0_truth==1)' % sel, wt))
    return res


def ApplyPhotonSplitting(hlist, components=['R', 'J', 'E']):
    res = []
    for label, sel, wt in hlist:
        if 'R' in components:
            res.append(('%s_R' % label, '(%s) && $p0_isprompt' % sel, wt))
        if 'J' in components:
            res.append(('%s_J' % label, '(%s) && $p0_isjet' % sel, wt))
        if 'E' in components:
            res.append(('%s_E' % label, '(%s) && $p0_iselec' % sel, wt))
    return res


def ApplySystWeightSplitting(hlist, wt_systs=[]):
    res = []
    for label, sel, wt in hlist:
        # res.append((label, sel, wt))
        for wt_expr, wt_name in wt_systs:
            res.append(('%s_%sUp' % (label, wt_name), sel, '(%s) * %s_hi' % (wt, wt_expr)))
            res.append(('%s_%sDown' % (label, wt_name), sel, '(%s) * %s_lo' % (wt, wt_expr)))
    return res


def ApplyScaleWeightSplitting(hlist):
    res = []
    for label, sel, wt in hlist:
        for i in xrange(6):
            res.append(('%s_scale_%i' % (label, i), sel, '(%s) * wt_sc_%i' % (wt, i)))
    return res


def ApplyPostFix(hlist, postfix):
    res = []
    for label, sel, wt in hlist:
        res.append(('%s_%s' % (label, postfix), sel, wt))
    return res


def StandardHists(node, var_list, binning, sel, wt, chn, manager, wt_systs=[], doFakes=True, doSysts=True):
    # Standard treatment for most MC samples
    for grp in groups:
        for P in groups[grp]:
            if grp == 'DY':
                hlist = ApplyPhotonSplitting(ApplyDYSplitting([(P, sel, wt)], components=['XZG']), components=['R'])
                hlist.extend(ApplyPhotonSplitting([(P, sel, wt)], components=['J', 'E']))
            elif grp == 'ZG':
                hlist = ApplyPhotonSplitting(ApplyDYSplitting([(P, sel, wt)], components=['IZG']), components=['R'])
            elif grp == 'TT':
                hlist = ApplyPhotonSplitting(ApplyTTSplitting([(P, sel, wt)], components=['XTTG']), components=['R'])
                hlist.extend(ApplyPhotonSplitting([(P, sel, wt)], components=['J', 'E']))
            elif grp == 'TTG':
                hlist = ApplyPhotonSplitting(ApplyTTSplitting([(P, sel, wt)], components=['ITTG']), components=['R'])
            else:
                hlist = ApplyPhotonSplitting([(P, sel, wt)])
            hlist_fw = []
            if doFakes:
                hlist_fw = ApplyPostFix(hlist, 'fw')
            hlist_init = list(hlist)
            if doSysts:
                hlist.extend(ApplySystWeightSplitting(hlist_init, wt_systs))
                if grp == 'WG':
                    hlist.extend(ApplyScaleWeightSplitting(hlist_init))
            # pprint(hlist)
            # pprint(hlist_fw)
            for H in hlist:
                node[H[0]] = Hist('TH1F', sample=P, var=var_list, binning=binning, sel=X.get(H[1]), wt=X.get(H[2]))
            for H in hlist_fw:
                node[H[0]] = Hist('TH1F', sample=P, var=var_list, binning=binning, sel=X.get(H[1], override={"sig_t": "$sig_l"}), wt=X.get('(%s) * wt_p0_fake' % H[2]))

    node['data_obs'] = Hist('TH1F', sample='data_obs_%s' % chn, var=var_list, binning=binning, sel=X.get(sel), wt='1')
    if doFakes:
        fake_sel = X.get(sel, override={"sig_t": "$sig_l"})
        node['data_fakes'] = Hist('TH1F', sample='data_obs_%s' % chn, var=var_list, binning=binning, sel=fake_sel, wt='wt_p0_fake')
        loose_fake_sel = X.get(sel, override={"photon_sel": "!p0_loose"})
        node['data_fakes_highpt'] = Hist('TH1F', sample='data_obs_%s' % chn, var=var_list, binning=binning, sel=loose_fake_sel, wt='wt_p0_highpt_fake')
        for ifake in xrange(1, 13):
            node['data_fakes_WeightStatSystBin%iUp' % ifake] = Hist('TH1F', sample='data_obs_%s' % chn, var=var_list, binning=binning, sel=fake_sel, wt='wt_p0_fake * (1. + wt_p0_fake_err * (wt_p0_fake_bin == %i))' % ifake)
            node['data_fakes_WeightStatSystBin%iDown' % ifake] = Hist('TH1F', sample='data_obs_%s' % chn, var=var_list, binning=binning, sel=fake_sel, wt='wt_p0_fake * (1. - wt_p0_fake_err * (wt_p0_fake_bin == %i))' % ifake)

def HistSum(label_list):
    return sum([node[X] for X in label_list[1:]], node[label_list[0]])


def CapNegativeBins(h):
    for ibin in xrange(1, h.GetNbinsX() + 1):
        if h.GetBinContent(ibin) < 0.:
            h.SetBinContent(ibin, 0.)


def MakeScaleEnvelope(node, label):
    h = node[label]
    h_hi = node[label].Clone()
    h_lo = node[label].Clone()
    h_alts = []
    for i in range(6):
        h_alts.append(node[label + '_scale_%i' % i])
    for ix in xrange(1, h.GetNbinsX() + 1):
        max_dev = max([abs(x.GetBinContent(ix) - h.GetBinContent(ix)) for x in h_alts])
        h_hi.SetBinContent(ix, h_hi.GetBinContent(ix) + max_dev)
        h_lo.SetBinContent(ix, max(0., h_lo.GetBinContent(ix) - max_dev))
    return h_hi, h_lo


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
X['iso_l'] = 'p0_chiso > 4 && p0_chiso < 10'
# X['iso_l'] = 'p0_chiso > 4 && p0_chiso < 10'
X['sig_t'] = '($p0_eb && p0_sigma < 0.01015) || ($p0_ee && p0_sigma < 0.0272)'
# X['sig_l'] = '($p0_eb && p0_sigma > 0.01015) || ($p0_ee && p0_sigma > 0.0272)'
X['sig_l'] = '($p0_eb && p0_sigma > 0.01100) || ($p0_ee && p0_sigma > 0.0300)'
X['photon_sel'] = 'p0_medium_noch && $iso_t && $sig_t'
# Selections for vetoing e->photon fakes
X['efake_veto_e'] = '!p0_haspix && p0_eveto && n_veto_e == 1 && n_veto_m == 0'
X['efake_veto_m'] = '!p0_haspix && p0_eveto && n_veto_m == 1 && n_veto_e == 0'

# Selection to veto (or select) mZ region
X['mZ_veto_m'] = '(l0p0_M < 70 || l0p0_M > 100)'
X['mZ_veto_e'] = '(l0p0_M < 70 || l0p0_M > 110)'
X['mZ_veto_inv_m'] = '(l0p0_M >= 70 && l0p0_M <= 100)'
X['mZ_veto_inv_e'] = '(l0p0_M >= 75 && l0p0_M <= 105)'

# Veto specific boxes in photon phi-eta where we have e->gamma spikes
if args.year == '2016':
    X['p0_phi_veto'] = 'p0_phi > 2.268 && p0_phi < 2.772 && p0_eta > -2.04 && p0_eta < 0.84'
if args.year == '2017':
    X['p0_phi_veto'] = 'p0_phi > 2.646 && p0_phi < 3.15 && p0_eta > -0.60 && p0_eta < 1.68'
if args.year == '2018':
    X['p0_phi_veto'] = 'p0_phi > 0.504 && p0_phi < 0.882 && p0_eta > -1.44 && p0_eta < 1.44'

# Analysis selection levels:
X['baseline_m_nopix'] ='metfilters==0 && l0_pdgid == 13 && l0_trg && n_pre_m>=1 && n_pre_p==1 && $photon_sel && l0_pt>30 && abs(l0_eta) < 2.4 && p0_pt>30 && abs(p0_eta) < 2.5 && l0p0_dr>0.7'
X['baseline_e_nopix'] ='metfilters==0 && l0_pdgid == 11 && l0_trg && n_pre_e>=1 && n_pre_p==1 && $photon_sel && l0_pt>35 && abs(l0_eta) < 2.5 && p0_pt>30 && abs(p0_eta) < 2.5 && l0p0_dr>0.7 && !$p0_phi_veto'
X['baseline_m'] ='$baseline_m_nopix && $efake_veto_m'
X['baseline_e'] ='$baseline_e_nopix && $efake_veto_e'
X['baseline_m_met'] ='$baseline_m && puppi_met>40'
X['baseline_e_met'] ='$baseline_e && puppi_met>40'
X['baseline_m_nomet'] ='$baseline_m && $mZ_veto_m'
X['baseline_e_nomet'] ='$baseline_e && $mZ_veto_e'
X['baseline_m_mZ_veto'] ='$baseline_m && $mZ_veto_m && puppi_met>40'
X['baseline_e_mZ_veto'] ='$baseline_e && $mZ_veto_e && puppi_met>40'

# Control regions
X['cr_Zmm'] ='l0_pdgid == 13 && l0_trg && n_pre_m==2 && l0_pt>30 && l1_pt>30 && met>0 && l0l1_os && l0l1_dr > 0.3'
X['cr_Zee'] ='l0_pdgid == 11 && l0_trg && n_pre_e==2 && l0_pt>30 && l1_pt>30 && met>0 && l0l1_os && l0l1_dr > 0.3'

# Common fiducial region:
if args.task == 'eft_region':
    # X['fid_m'] ='l0_pdgid == 13 && l0_trg && n_pre_m==1 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>80 && met>80 && p0_pt>150 && l0p0_dr>3.0 && $efake_veto_m'
    # X['fid_e'] ='l0_pdgid == 11 && l0_trg && n_pre_e==1 && n_pre_p==1 && p0_medium_noch && $iso_t && $sig_t && l0_pt>80 && met>80 && p0_pt>150 && l0p0_dr>3.0 && $mZ_veto_e && $efake_veto_e'
    X['fid_m'] ='$baseline_m_mZ_veto && l0_pt>80 && puppi_met>80 && p0_pt>150 && l0p0_dr>3.0'
    X['fid_e'] ='$baseline_e_mZ_veto && l0_pt>80 && puppi_met>80 && p0_pt>150 && l0p0_dr>3.0'
    X['gen_core'] = 'gen_l0_pt>80 && abs(gen_l0_eta)<2.5 && gen_met>40 && gen_p0_pt>150 && abs(gen_p0_eta)<2.5 && gen_l0p0_dr>3.0 && lhe_frixione'
    X['gen_acc'] = '$gen_core && gen_met>80'
    X['gen_acc_met1'] = '$gen_core && gen_met>40 && gen_met<=80'
    X['gen_ooa'] = '!$gen_core'

if args.task == 'fid_region':
    X['fid_m'] ='$baseline_m_mZ_veto'
    X['fid_e'] ='$baseline_e_mZ_veto'
    X['gen_core'] = 'gen_l0_pt>30 && abs(gen_l0_eta)<2.5 && gen_met>0 && gen_p0_pt>30 && abs(gen_p0_eta)<2.5 && gen_l0p0_dr>0.7 && lhe_frixione'
    X['gen_acc'] = '$gen_core && gen_met>40'
    X['gen_acc_met1'] ='$gen_core && gen_met>0 && gen_met<=40'
    X['gen_ooa'] = '!$gen_core'

# Event weights
X['baseline_wt'] = 'wt_def*wt_pu*wt_l0*wt_trg_l0*wt_p0*wt_p0_e_fake*wt_pf'
X['cr_zll_wt'] = 'wt_def*wt_pu*wt_l0*wt_trg_l0*wt_l1*wt_pf'


if args.task == 'eft_region' or args.task == 'fid_region':
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

    do_cats['e'].extend(['fid_e'])
    do_cats['m'].extend(['fid_m'])

    X['p_gen_ooa'] = 'gen_l0_q==+1 && $gen_ooa'
    X['n_gen_ooa'] = 'gen_l0_q==-1 && $gen_ooa'
    X['x_gen_ooa'] = '$gen_ooa'

    X['p_gen_acc'] = 'gen_l0_q==+1 && $gen_acc'
    X['n_gen_acc'] = 'gen_l0_q==-1 && $gen_acc'
    X['x_gen_acc'] = '$gen_acc'

    X['p_gen_acc_met1'] = 'gen_l0_q==+1 && $gen_acc_met1'
    X['n_gen_acc_met1'] = 'gen_l0_q==-1 && $gen_acc_met1'
    X['x_gen_acc_met1'] = '$gen_acc_met1'

    for i in range(len(pt_bins_min)):
        # Reconstructed categories
        X['x_e_%i' % i] = '$fid_e && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['x_m_%i' % i] = '$fid_m && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])

        X['p_e_%i' % i] = '$fid_e && l0_q==+1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_e_%i' % i] = '$fid_e && l0_q==-1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['p_m_%i' % i] = '$fid_m && l0_q==+1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_m_%i' % i] = '$fid_m && l0_q==-1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        if args.task == 'eft_region':
            do_cats['e'].extend(['p_e_%i' % i, 'n_e_%i' % i])
            do_cats['m'].extend(['p_m_%i' % i, 'n_m_%i' % i])
        if args.task == 'fid_region':
            do_cats['e'].extend(['x_e_%i' % i])
            do_cats['m'].extend(['x_m_%i' % i])

        # Gen level selections
        X['x_gen_%i' % (i)] = '$gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['x_gen_met1_%i' % (i)] = '$gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['p_gen_%i' % (i)] = '$p_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_gen_%i' % (i)] = '$n_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['p_gen_met1_%i' % (i)] = '$p_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])
        X['n_gen_met1_%i' % (i)] = '$n_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i])

        for j in range(len(phi_bins_min)):
            X['p_gen_%i_%i' % (i, j)] = '$p_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['n_gen_%i_%i' % (i, j)] = '$n_gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['x_gen_%i_%i' % (i, j)] = '$gen_acc && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['p_gen_met1_%i_%i' % (i, j)] = '$p_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['n_gen_met1_%i_%i' % (i, j)] = '$n_gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])
            X['x_gen_met1_%i_%i' % (i, j)] = '$gen_acc_met1 && gen_p0_pt>=%f && gen_p0_pt<%f && %s >= %f && %s < %f' % (pt_bins_min[i], pt_bins_max[i], phi_var, phi_bins_min[j], phi_var, phi_bins_max[j])

    drawvars = [
        # ('gen_p0_pt', BinningFromStr(eft_defaults['pt_bins'])),
        ('p0_pt', BinningFromStr(eft_defaults['pt_bins'])),
        # (eft_defaults['phi_var'], BinningFromStr(eft_defaults['phi_bins'])),
        (eft_defaults['phi_var_obs'], BinningFromStr(eft_defaults['phi_bins_obs']))
    ]

    for chn in ['e', 'm']:
        for sel in do_cats[chn]:
            for var, binning in drawvars:
                xnode = hists[chn][sel][var]
                StandardHists(hists[chn][sel][var], var_list=[var], binning=binning, sel=('$' + sel), wt='$baseline_wt', chn=chn, manager=X, wt_systs=main_wt_systs[chn], doSysts=(args.syst is None))
                wg_scale_hlist = []
                wg_syst_hlist = []
                wg_hlist = []
                for chg in ['p', 'n', 'x']:
                    hlist = [
                        ('WG_ooa_%s' % chg, '$' + sel + ' && $%s_gen_ooa' % chg, '$baseline_wt'),
                        ('WG_main_%s' % chg, '$' + sel + ' && $%s_gen_acc' % chg, '$baseline_wt'),
                        ('WG_met1_%s' % chg, '$' + sel + ' && $%s_gen_acc_met1' % chg, '$baseline_wt'),
                    ]
                    wg_hlist.extend(hlist)
                    wg_scale_hlist.extend(hlist)
                    wg_syst_hlist.extend(hlist)
                    # xnode['WG_ooa_%s' % chg] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_ooa' % chg), wt=X.get('$baseline_wt'))
                    # xnode['WG_main_%s' % chg] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_acc' % chg), wt=X.get('$baseline_wt'))
                    # xnode['WG_met1_%s' % chg] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_acc_met1' % chg), wt=X.get('$baseline_wt'))
                if sel[-1].isdigit():
                    chg = sel[0]
                    for i in range(len(pt_bins_min)):
                        hlist = [
                            ('WG_main_%s_%i' % (chg, i), '$' + sel + ' && $%s_gen_%i' % (chg, i), '$baseline_wt'),
                            ('WG_met1_%s_%i' % (chg, i), '$' + sel + ' && $%s_gen_met1_%i' % (chg, i), '$baseline_wt'),
                        ]
                        wg_hlist.extend(hlist)
                        wg_syst_hlist.extend(hlist)

                        # xnode['WG_main_%s_%i' % (chg, i)] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_%i' % (chg, i)), wt=X.get('$baseline_wt'))
                        # xnode['WG_met1_%s_%i' % (chg, i)] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_met1_%i' % (chg, i)), wt=X.get('$baseline_wt'))
                    this_pt_bin = int(sel[-1])
                    do_pt_bins_in_phi = [this_pt_bin]
                    if this_pt_bin > 0:
                        do_pt_bins_in_phi.append(this_pt_bin - 1)
                    if this_pt_bin < len(pt_bins_min) - 1:
                        do_pt_bins_in_phi.append(this_pt_bin + 1)
                    for i in do_pt_bins_in_phi:
                        for j in range(len(phi_bins_min)):
                            hlist = [
                                ('WG_main_%s_%i_%i' % (chg, i, j), '$' + sel + ' && $%s_gen_%i_%i' % (chg, i, j), '$baseline_wt'),
                                ('WG_met1_%s_%i_%i' % (chg, i, j), '$' + sel + ' && $%s_gen_met1_%i_%i' % (chg, i, j), '$baseline_wt'),
                            ]
                            wg_hlist.extend(hlist)
                            wg_scale_hlist.extend(hlist)
                            wg_syst_hlist.extend(hlist)
                            # xnode['WG_main_%s_%i_%i' % (chg, i, j)] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_%i_%i' % (chg, i, j)), wt=X.get('$baseline_wt'))
                            # xnode['WG_met1_%s_%i_%i' % (chg, i, j)] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('$' + sel + ' && $%s_gen_met1_%i_%i' % (chg, i, j)), wt=X.get('$baseline_wt'))
                if args.syst is None:
                    wg_hlist.extend(ApplySystWeightSplitting(wg_syst_hlist, main_wt_systs[chn]))
                    wg_hlist.extend(ApplyScaleWeightSplitting(wg_scale_hlist))
                for H in wg_hlist:
                    xnode[H[0]] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get(H[1]), wt=X.get(H[2]))

        # Rest of this only for the nominal run
        if args.syst is not None:
            continue
        pdgid = '11' if chn == 'e' else '13'
        for var, binning in drawvars:
            hists[chn]['XS'][var]['XS_WG_p_%s_acc' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $p_gen_acc' % pdgid), wt=X.get('wt_def'))
            hists[chn]['XS'][var]['XS_WG_n_%s_acc' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $n_gen_acc' % pdgid), wt=X.get('wt_def'))
            hists[chn]['XS'][var]['XS_WG_x_%s_acc' % chn] = Hist('TH1F', sample=wg_sample, var=[var], binning=binning, sel=X.get('gen_pdgid==%s && $x_gen_acc' % pdgid), wt=X.get('wt_def'))
        hists[chn]['XS']['2D']['XS_WG_p_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $p_gen_acc' % pdgid), wt=X.get('wt_def'))
        hists[chn]['XS']['2D']['XS_WG_n_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $n_gen_acc' % pdgid), wt=X.get('wt_def'))
        hists[chn]['XS']['2D']['XS_WG_x_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $x_gen_acc' % pdgid), wt=X.get('wt_def'))
        # Also make a copy that does get normalised to the lumi, for doing
        hists[chn]['XS']['2D']['WG_p_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $p_gen_acc' % pdgid), wt=X.get('wt_def'))
        hists[chn]['XS']['2D']['WG_n_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $n_gen_acc' % pdgid), wt=X.get('wt_def'))
        hists[chn]['XS']['2D']['WG_x_%s_acc' % chn] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $x_gen_acc' % pdgid), wt=X.get('wt_def'))
        for i_sc in range(6):
            hists[chn]['XS']['2D']['XS_WG_p_%s_acc__sc_%i' % (chn, i_sc)] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $p_gen_acc' % pdgid), wt=X.get('wt_def * wt_sc_%i' % i_sc))
            hists[chn]['XS']['2D']['XS_WG_n_%s_acc__sc_%i' % (chn, i_sc)] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $n_gen_acc' % pdgid), wt=X.get('wt_def * wt_sc_%i' % i_sc))
            hists[chn]['XS']['2D']['XS_WG_x_%s_acc__sc_%i' % (chn, i_sc)] = Hist('TH2F', sample=wg_sample, var=['gen_p0_pt', phi_var], binning=(BinningFromStr(eft_defaults['pt_bins']) + BinningFromStr(eft_defaults['phi_bins'])), sel=X.get('gen_pdgid==%s && $x_gen_acc' % pdgid), wt=X.get('wt_def * wt_sc_%i' % i_sc))

if args.task in ['baseline', 'electron_fakes']:
    wt_systs = dict(main_wt_systs)
    if args.task == 'electron_fakes':
        pt_bins_min = [30, 35, 40, 60, 100]
        pt_bins_max = [35, 40, 60, 100, 200]
        eta_bins_min = [0, 1.4442]
        eta_bins_max = [1.4442, 2.5]
        X['baseline_wt'] = 'wt_def*wt_pu*wt_l0*wt_trg_l0*wt_p0*wt_pf'
        wt_systs = {'e': [], 'm': []}

        for e_min, e_max in zip(eta_bins_min, eta_bins_max):
            for p_min, p_max in zip(pt_bins_min, pt_bins_max):
                label = '%g_%g_%g_%g' % (e_min, e_max, p_min, p_max)
                label = label.replace('.', 'p')
                X['fail_%s' % label] ='$baseline_e_nopix && $mZ_veto_inv_e && !$efake_veto_e && abs(p0_eta) >= %g && abs(p0_eta) < %g && p0_pt >= %g && p0_pt < %g' % (e_min, e_max, p_min, p_max)
                X['pass_%s' % label] ='$baseline_e_nopix && $mZ_veto_inv_e && $efake_veto_e && abs(p0_eta) >= %g && abs(p0_eta) < %g && p0_pt >= %g && p0_pt < %g' % (e_min, e_max, p_min, p_max)
                X['tot_%s' % label] ='$baseline_e_nopix && $mZ_veto_inv_e && abs(p0_eta) >= %g && abs(p0_eta) < %g && p0_pt >= %g && p0_pt < %g' % (e_min, e_max, p_min, p_max)
                do_cats['e'].append('fail_%s' % label)
                do_cats['e'].append('pass_%s' % label)
                do_cats['e'].append('tot_%s' % label)

    do_cats['e'].extend(['baseline_e_nopix', 'baseline_e', 'baseline_e_mZ_veto', 'baseline_e_met', 'baseline_e_nomet', 'cr_Zee'])
    do_cats['m'].extend(['baseline_m_nopix', 'baseline_m', 'baseline_m_mZ_veto', 'baseline_m_met', 'baseline_m_nomet', 'cr_Zmm'])

    cat_to_wt = {
        'cr_Zee': '$cr_zll_wt',
        'cr_Zmm': '$cr_zll_wt'
    }

    drawvars = [
        ('n_vtx', (30, 0., 60.)),
        ('n_veto_m', (3, -0.5, 2.5)),
        ('n_veto_e', (3, -0.5, 2.5)),
        ('l0_pt', (30, 0., 150.)),
        ('l0_eta', (20, -3.0, 3.0)),
        ('l0_phi', (20, -3.15, 3.15)),
        # ('l0_iso', (40, 0, 2.0)),
        ('l1_pt', (40, 0., 150.)),
        ('l0l1_M', (60, 60, 120)),
        ('l0l1_pt', (80, 0, 200)),
        # ('l0l1_dr', (20, 0., 5.)),
        ('met', (20, 0., 200.)),
        # ('tk_met', (40, 0., 200.)),
        # ('tk_met_phi', (20, -3.15, 3.15)),
        ('puppi_met', (20, 0., 200.)),
        ('p0_pt', [0, 10, 20, 30, 40, 50, 60, 80, 100, 120, 160, 200, 250, 300, 400, 500, 600, 850, 1200]),
        ('p0_eta', (20, -3.0, 3.0)),
        ('p0_phi', (30, -3.15, 3.15)),
        ('l0p0_dr', (20, 0., 5.)),
        ('l0p0_M', (40, 0, 200)),
        ('p0_chiso', (40, 0, 2.0)),
        # ('p0_neiso', (40, 0, 20.0)),
        # ('p0_phiso', (40, 0, 20.0)),
        # ('p0_hovere', (20, 0., 0.5)),
        ('p0_sigma', (60, 0., 0.06)),
        ('p0_haspix', (2, -0.5, 1.5)),
        ('p0_eveto', (2, -0.5, 1.5)),
        ('p0_truth', (7, -0.5, 6.5)),
        ('wt_def', (100, 0.8, 1.2)),
        ('wt_pu', (100, 0.5, 1.5)),
        ('wt_l0', (100, 0.8, 1.2)),
        ('wt_l1', (100, 0.8, 1.2)),
        ('wt_trg_l0', (100, 0.8, 1.2)),
        ('wt_p0', (100, 0.8, 1.2)),
        ('wt_p0_e_fake', (100, 0, 4))
    ]

    if args.syst == None:
        drawvars.extend([
            ('l0met_mt', (30, 0., 200.)),
            ('l1_eta', (20, -3.0, 3.0)),
            ('met_phi', (20, -3.15, 3.15)),
            ('puppi_met_phi', (20, -3.15, 3.15))
        ])


    if args.task == 'electron_fakes':
        drawvars = [
            ('l0p0_M', (40, 60, 120)),
        ]

    for chn in ['e', 'm']:
        for sel in do_cats[chn]:
            for var, binning in drawvars:
                wt = '$baseline_wt'
                if sel in cat_to_wt:
                    wt = cat_to_wt[sel]
                if var.startswith('wt_'):
                    wt = '1.0'
                if 'mZ_veto' in sel:
                    actual_wt_systs = wt_systs[chn]
                else:
                    actual_wt_systs = list()
                StandardHists(hists[chn][sel][var], var_list=[var], binning=binning, sel=('$' + sel), wt=wt, chn=chn, manager=X, wt_systs=actual_wt_systs, doSysts=(args.syst is None))
            for sample in samples:
                hists['2D'][chn][sel]['l0_eta_phi'][sample] = Hist('TH2F', sample=sample, var=['l0_eta', 'l0_phi'], binning=(50, -3, 3, 50, -3.15, 3.15), sel=X.get('$' + sel), wt=X.get(wt))
                hists['2D'][chn][sel]['p0_eta_phi'][sample] = Hist('TH2F', sample=sample, var=['p0_eta', 'p0_phi'], binning=(50, -3, 3, 50, -3.15, 3.15), sel=X.get('$' + sel), wt=X.get(wt))


if args.task in ['photon_fakes', 'photon_fakes2']:
    if args.task == 'photon_fakes2':
        X['baseline_m_nopix'] ='l0_pdgid == 13 && l0_trg && n_pre_m==1 && n_pre_p==1 && $photon_sel && l0_pt>30 && abs(l0_eta) < 2.4 && p0_pt>100 && abs(p0_eta) < 2.5 && l0p0_dr>0.7'
        X['baseline_e_nopix'] ='l0_pdgid == 11 && l0_trg && n_pre_e==1 && n_pre_p==1 && $photon_sel && l0_pt>35 && abs(l0_eta) < 2.5 && p0_pt>100 && abs(p0_eta) < 2.5 && l0p0_dr>0.7 && !$p0_phi_veto'
        X['baseline_m_mZ_veto'] ='$baseline_m && $mZ_veto_m && puppi_met<40'
        X['baseline_e_mZ_veto'] ='$baseline_e && $mZ_veto_e && puppi_met<40'
        extra_regions = [
            ('_pt_100_150', 'p0_pt>100 && p0_pt<=150'),
            ('_pt_150_200', 'p0_pt>150 && p0_pt<=200'),
            ('_pt_200_300', 'p0_pt>200 && p0_pt<=300'),
            ('_pt_300_400', 'p0_pt>300 && p0_pt<=400'),
            ('_pt_400_600', 'p0_pt>400 && p0_pt<=600'),
        ]
    if args.task == 'photon_fakes':
        extra_regions = [
            ('_pt_30_40', 'p0_pt>30 && p0_pt<=40'),
            ('_pt_40_60', 'p0_pt>40 && p0_pt<=60'),
            ('_pt_60_100', 'p0_pt>60 && p0_pt<=100'),
            ('_pt_100_200', 'p0_pt>100 && p0_pt<=200'),
        ]


    X['baseline_m_mZ_veto'] = X.get('$baseline_m_mZ_veto', override={"sig_t": "1", "iso_t": "1"})
    X['baseline_e_mZ_veto'] = X.get('$baseline_e_mZ_veto', override={"sig_t": "1", "iso_t": "1"})
    X['barrel_m'] = '$baseline_m_mZ_veto && $p0_eb'
    X['endcap_m'] = '$baseline_m_mZ_veto && $p0_ee'
    X['barrel_e'] = '$baseline_e_mZ_veto && $p0_eb'
    X['endcap_e'] = '$baseline_e_mZ_veto && $p0_ee'

    for chn in ['e', 'm']:
        for S in ['barrel_%s' % chn, 'endcap_%s' % chn]:
            for POST, EXTRA in [
              ('', '1'),
            ] + extra_regions:
                X['%s_iso_l%s' % (S, POST)] = '$%s && $iso_l && (%s)' % (S, EXTRA)
                X['%s_sig_l%s' % (S, POST)] = '$%s && $sig_l && (%s)' % (S, EXTRA)
                X['%s_iso_t%s' % (S, POST)] = '$%s && $iso_t && (%s)' % (S, EXTRA)
                X['%s_sig_t%s' % (S, POST)] = '$%s && $sig_t && (%s)' % (S, EXTRA)
                X['%s_iso_l_sig_t%s' % (S, POST)] = '$%s && $iso_l && $sig_t && (%s)' % (S, EXTRA)
                X['%s_iso_l_sig_l%s' % (S, POST)] = '$%s && $iso_l && $sig_l && (%s)' % (S, EXTRA)
                X['%s_iso_t_sig_l%s' % (S, POST)] = '$%s && $iso_t && $sig_l && (%s)' % (S, EXTRA)
                X['%s_iso_t_sig_t%s' % (S, POST)] = '$%s && $iso_t && $sig_t && (%s)' % (S, EXTRA)
                do_cats[chn].extend(
                    ['%s_iso_l%s' % (S, POST),
                     '%s_sig_l%s' % (S, POST),
                     '%s_iso_t%s' % (S, POST),
                     '%s_sig_t%s' % (S, POST),
                     '%s_iso_l_sig_t%s' % (S, POST),
                     '%s_iso_l_sig_l%s' % (S, POST),
                     '%s_iso_t_sig_l%s' % (S, POST),
                     '%s_iso_t_sig_t%s' % (S, POST)])

    drawvars = [
        ('l0met_mt', (30, 0., 200.)),
        ('l0_pt', (40, 0., 150.)),
        ('l0_eta', (20, -3.0, 3.0)),
        # ('l0l1_M', (40, 60, 120)),
        ('met', (20, 0., 200.)),
        ('p0_pt', (120, 0, 600.)),
        ('p0_eta', (20, -3.0, 3.0)),
        # ('l0p0_dr', (20, 0., 5.)),
        # ('l0p0_M', (20, 60, 120)),
        ('p0_chiso', (40, 0, 20.0)),
        # ('p0_neiso', (40, 0, 20.0)),
        # ('p0_phiso', (40, 0, 20.0)),
        # ('p0_hovere', (20, 0., 0.5)),
        ('p0_sigma', (60, 0., 0.06)),
        ('p0_haspix', (2, -0.5, 1.5)),
        ('p0_truth', (7, -0.5, 6.5)),
    ]

    for chn in ['e', 'm']:
        for sel in do_cats[chn]:
            for var, binning in drawvars:
                doFakes = ('iso_t_sig_t' in sel)
                StandardHists(hists[chn][sel][var], var_list=[var], binning=binning, sel=('$' + sel), wt='$baseline_wt', chn=chn, manager=X, wt_systs=[], doFakes=doFakes, doSysts=False)

MultiDraw(hists, samples, tname, mt_cores=4, mt_thresh=1)

with open('input/cfg_wgamma_%s_v4.json' % year) as jsonfile:
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

# fout = ROOT.TFile('output_%s_%s_%s.root' % (year, args.task, args.label), 'RECREATE')
# NodeToTDir(fout, hists, args.syst)
# fout.Close()
# sys.exit(0)

# fin = ROOT.TFile('output_%s_%s_%s.root' % (year, args.task, args.label))
# TDirToNode(fin, node=hists)
doCleanup = False

allProcs = []
for sa in groups['WG']:
    allProcs.append('%s_R' % sa)
for sa in groups['DY']:
    allProcs.append('%s_XZG_R' % sa)
    allProcs.append('%s_E' % sa)
for sa in groups['ZG']:
    allProcs.append('%s_IZG_R' % sa)
for sa in groups['TT']:
    allProcs.append('%s_XTTG_R' % sa)
    allProcs.append('%s_E' % sa)
for sa in groups['TTG']:
    allProcs.append('%s_ITTG_R' % sa)
for sa in groups['VV']:
    allProcs.append('%s_R' % sa)
    allProcs.append('%s_E' % sa)
for sa in groups['GG']:
    allProcs.append('%s_R' % sa)
    allProcs.append('%s_E' % sa)
allProcs.append('data_fakes_highpt')

for path, node in hists.ListNodes(withObjects=True):
    print path
    if path.startswith('2D'):
        continue
    if 'XS' in path.split('/'):
        continue

    # CheckBinErrors(node, allProcs)

    for grp, grplist in groups.iteritems():
        suffixes = [g.replace(grplist[0] + '_', '') for g in node.d.keys() if g.startswith(grplist[0] + '_')]
        scale_variations = []
        for suf in suffixes:
            # print grp, grplist, suf
            if len(grplist) > 1:
                sum_list = [g + '_' + suf for g in grplist]
                node[grp + '_' + suf] = HistSum(sum_list)
                if doCleanup:
                    for hname in sum_list:
                        del node.d[hname]
            if '_scale_0' in suf:
                scale_hi, scale_lo = MakeScaleEnvelope(node, grp + '_' + suf.replace('_scale_0', ''))
                node[grp + '_' + suf.replace('_scale_0', '_QCDScaleUp')] = scale_hi
                node[grp + '_' + suf.replace('_scale_0', '_QCDScaleDown')] = scale_lo
                # scale_variations.append(suf.replace('_scale_0', ''))



    node['Total_R'] = HistSum(['WG_R', 'TT_XTTG_R', 'TTG_ITTG_R', 'DY_XZG_R', 'ZG_IZG_R', 'VV_R'])
    node['Total_E'] = HistSum(['W_E', 'TT_E', 'DY_E', 'VV_E'])
    node['Total_J'] = HistSum(['W_J', 'TT_J', 'DY_J', 'VV_J'])

    if 'data_fakes' in node.d:
        node['Total_R_fw'] = HistSum(['WG_R_fw', 'TT_XTTG_R_fw', 'TTG_ITTG_R_fw', 'DY_XZG_R_fw', 'ZG_IZG_R_fw', 'VV_R_fw'])
        node['Total_E_fw'] = HistSum(['W_E_fw', 'TT_E_fw', 'DY_E_fw', 'VV_E_fw'])
        node['data_fakes_sub'] = node['data_fakes'] - (node['Total_R_fw'] + node['Total_E_fw'])
        CapNegativeBins(node['data_fakes_sub'])
        for ifake in xrange(1, 13):
            node['data_fakes_sub_WeightStatSystBin%iUp' % ifake] = node['data_fakes_WeightStatSystBin%iUp' % ifake] - (node['Total_R_fw'] + node['Total_E_fw'])
            node['data_fakes_sub_WeightStatSystBin%iDown' % ifake] = node['data_fakes_WeightStatSystBin%iDown' % ifake] - (node['Total_R_fw'] + node['Total_E_fw'])
            CapNegativeBins(node['data_fakes_sub_WeightStatSystBin%iUp' % ifake])
            CapNegativeBins(node['data_fakes_sub_WeightStatSystBin%iDown' % ifake])

fout = ROOT.TFile('output_%s_%s_%s.root' % (year, args.task, args.label), 'RECREATE')
NodeToTDir(fout, hists, args.syst)
fout.Close()
