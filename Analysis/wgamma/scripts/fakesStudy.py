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
parser.add_argument('--year', default='2016', choices=['2016', '2017', '2018'])
parser.add_argument('--indir', default='output/130818/wgamma_2016_v2/WGamma/')

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

year = args.year

remaps = {
    "2016": {
        'W': 'WJetsToLNu-stitched',
    },
    "2017": {

    },
    "2018": {

    }
}

remap = remaps[year]

# ttgammas = ['TTG_DL', 'TTG_Had', 'TTG_SL_T', 'TTG_SL_Tbar']
# zg_samples = ['ZG']
# if year == '2016':
#     dibosons = ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
#     tt_samples = ['TT']
# if year == '2017':
#     dibosons = ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
#     tt_samples = ['TT_SL', 'TT_Had', 'TT_DL']
# if year == '2018':
#     dibosons = ['WW', 'WZ', 'ZZ', 'ST_s', 'ST_t_antitop', 'ST_t_top', 'ST_tW_antitop', 'ST_tW_top']
    # tt_samples = ['TT_SL', 'TT_Had', 'TT_DL']

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')


def AddPhotonSplitting(node, label, sample, var_list, binning, sel, wt):
    node['%s_R' % label] = Hist('TH1D', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && $p0_isprompt' % sel), wt=X.get(wt))
    node['%s_F' % label] = Hist('TH1D', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && $p0_isfake' % sel), wt=X.get(wt))
    node['%s_J' % label] = Hist('TH1D', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && $p0_isjet' % sel), wt=X.get(wt))
    node['%s_E' % label] = Hist('TH1D', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && $p0_iselec' % sel), wt=X.get(wt))


def AddDYSplitting(node, label, sample, var_list, binning, sel, wt):
    node[label + '_XZG'] = Hist('TH1D', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && !(gen_is_zg && p0_truth==1)' % sel), wt=X.get(wt))
    node[label + '_IZG'] = Hist('TH1D', sample=sample, var=var_list, binning=binning, sel=X.get('(%s) && (gen_is_zg && p0_truth==1)' % sel), wt=X.get(wt))
    AddPhotonSplitting(node, label + '_XZG', sample, var_list, binning, '(%s) && !(gen_is_zg && p0_truth==1)' % sel, wt)
    AddPhotonSplitting(node, label + '_IZG', sample, var_list, binning, '(%s) && (gen_is_zg && p0_truth==1)' % sel, wt)

def HistSum(label_list):
    return sum([node[X] for X in label_list[1:]], node[label_list[0]])


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
X['iso_l'] = 'p0_chiso > 3 && p0_chiso < 10'
X['sig_t'] = '($p0_eb && p0_sigma < 0.01015) || ($p0_ee && p0_sigma < 0.0272)'
X['sig_l'] = '($p0_eb && p0_sigma > 0.01100) || ($p0_ee && p0_sigma > 0.0272)'

# Selections for vetoing e->photon fakes
X['efake_veto_e'] = '!p0_haspix && p0_eveto && n_veto_e == 1 && n_veto_m == 0'
X['efake_veto_m'] = '!p0_haspix && p0_eveto && n_veto_m == 1 && n_veto_e == 0'

# Analysis selection levels:
X['full_m'] ='l0_pdgid == 13 && n_pre_m==1 && l0_iso<0.15 && n_pre_p==1 && l0_pt>20 && met>20 && p0_pt>30 && l0p0_dr>0.7 && $p0_isjet && $efake_veto_m'
X['full_e'] ='l0_pdgid == 11 && n_pre_e==1 &&                n_pre_p==1 && l0_pt>20 && met>20 && p0_pt>30 && l0p0_dr>0.7 && $p0_isjet && $efake_veto_e'

X['baseline_m_nopix'] ='l0_pdgid == 13 && n_pre_m==1 && l0_iso<0.15 && n_pre_p==1 && p0_medium_noch && l0_pt>20 && met>20 && p0_pt>30 && l0p0_dr>0.7'
X['baseline_e_nopix'] ='l0_pdgid == 11 && n_pre_e==1                && n_pre_p==1 && p0_medium_noch && l0_pt>20 && met>20 && p0_pt>30 && l0p0_dr>0.7'
X['baseline_m'] ='$baseline_m_nopix && $efake_veto_m'
X['baseline_e'] ='$baseline_e_nopix && $efake_veto_e'
X['baseline_m_mZ_veto'] ='$baseline_m'
X['baseline_e_mZ_veto'] ='$baseline_e'

# Event weights
X['baseline_wt'] = 'wt_def*wt_pu*wt_l0*wt_trg_l0*wt_p0*wt_p0_e_fake*wt_pf'


p0_vars = ['p0_chiso', 'p0_neiso', 'p0_phiso', 'p0_hovere', 'p0_sigma']
p0_vars_binning = {
    'p0_chiso': (80, 0, 20.0),
    'p0_neiso': (80, 0, 20.0),
    'p0_phiso': (80, 0, 20.0),
    'p0_hovere': (80, 0, 0.15),
    'p0_sigma': (160, 0., 0.03),
}
p0_cuts = {
    'p0_chiso': 1.141,
    'p0_neiso': 1.7182,
    'p0_phiso': 2.220595,
    'p0_hovere': 0.02197,
    'p0_sigma': 0.01015
}

if args.task == 'photon_fakes':
    X['barrel_m'] = '$baseline_m_mZ_veto && $p0_eb'
    X['endcap_m'] = '$baseline_m_mZ_veto && $p0_ee'
    X['barrel_e'] = '$baseline_e_mZ_veto && $p0_eb'
    X['endcap_e'] = '$baseline_e_mZ_veto && $p0_ee'
    X['full_barrel_m'] = '$full_m && $p0_eb && p0_pt > 40 && p0_pt < 60'
    X['full_endcap_m'] = '$full_m && $p0_ee && p0_pt > 40 && p0_pt < 60'
    X['full_barrel_e'] = '$full_e && $p0_eb && p0_pt > 40 && p0_pt < 60'
    X['full_endcap_e'] = '$full_e && $p0_ee && p0_pt > 40 && p0_pt < 60'

    do_cats['e'].extend(['full_e', 'full_barrel_e', 'baseline_e_mZ_veto', 'barrel_e', 'endcap_e'])
    do_cats['m'].extend(['full_m', 'full_barrel_m', 'baseline_m_mZ_veto', 'barrel_m', 'endcap_m'])

    for chn in ['e', 'm']:
        for S in ['barrel_%s' % chn, 'endcap_%s' % chn]:
            X['%s_iso_l' % S] = '$%s && $iso_l' % S
            X['%s_sig_l' % S] = '$%s && $sig_l' % S
            X['%s_iso_l_sig_t' % S] = '$%s && $iso_l && $sig_t' % S
            X['%s_iso_l_sig_l' % S] = '$%s && $iso_l && $sig_l' % S
            X['%s_iso_t_sig_l' % S] = '$%s && $iso_t && $sig_l' % S
            X['%s_iso_t_sig_t' % S] = '$%s && $iso_t && $sig_t' % S
            X['%s_p0_medium' % S] = '$%s && p0_medium' % S
            do_cats[chn].extend(['%s_iso_l' % S, '%s_sig_l' % S, '%s_iso_l_sig_t' % S, '%s_iso_l_sig_l' % S, '%s_iso_t_sig_l' % S, '%s_iso_t_sig_t' % S, '%s_p0_medium' % S])

    drawvars = [
        # ('l0met_mt', (30, 0., 200.)),
        # ('l0_pt', (40, 0., 150.)),
        # ('l0_eta', (20, -3.0, 3.0)),
        # ('l0l1_M', (40, 60, 120)),
        # ('met', (20, 0., 200.)),
        ('p0_pt', [30, 40, 60, 100, 200, 300]),
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
                for P in ['W']:
                    AddPhotonSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')
                # for P in ['DY'] + zg_samples:
                #     AddDYSplitting(hists[chn][sel][var], P, P, [var], binning, '$' + sel, '$baseline_wt')

                #     hists[chn][sel][var]['data_obs'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning, sel=X.get('$' + sel), wt=X.get('$baseline_wt'))
                #     if 'iso_t_sig_t' in sel:
                #         hists[chn][sel][var]['data_fakes'] = Hist('TH1F', sample='data_obs_%s' % chn, var=[var], binning=binning,
                #             sel=X.get('$' + sel, override={"sig_t": "$sig_l"}, printlevel=0),
                #             wt=X.get('$baseline_wt * wt_p0_fake'))
            if not sel.startswith('full'):
                continue
            for ivar in xrange(0, len(p0_vars)):
                for jvar in xrange(ivar + 1, len(p0_vars)):
                    xname = p0_vars[ivar]
                    yname = p0_vars[jvar]
                    label = '%s_vs_%s' % (xname, yname)
                    cutstr = []
                    for zvar in p0_vars:
                        if zvar in [xname, yname]: continue
                        cutstr.append('%s < %f' % (zvar, p0_cuts[zvar]))
                    cutstr = ' && '.join(cutstr)
                    print xname, yname, cutstr
                    hists['2D'][chn][sel][label][sample] = Hist('TH2F',
                        sample=sample,
                        var=[xname, yname],
                        binning=(p0_vars_binning[xname] + p0_vars_binning[yname]),
                        sel=X.get('(%s) && $%s' %(cutstr, sel), printlevel=0),
                        wt=X.get('$baseline_wt')
                    )

preload = False

if not preload:
    MultiDraw(hists, samples, tname, mt_cores=4, mt_thresh=1)

    with open('input/cfg_wgamma_%s_v3.json' % year) as jsonfile:
        cfg = json.load(jsonfile)
        sample_cfg = cfg['samples']

    for sample in samples:
        f = ROOT.TFile.Open(samples[sample])
        sample_cfg[remap[sample]]['events'] = f.Get('counters').GetBinContent(2)
        f.Close()

    for chn in ['e', 'm', '2D']:
        for path, hname, obj in hists[chn].ListObjects():
            name = obj.sample
            if 'data_obs' not in name:
                tgt_lumi = 35916.764758
                # tgt_lumi = sample_cfg[remap['data_obs_%s' % chn]]['lumi']
                events = sample_cfg[remap[name]]['events']
                xsec = sample_cfg[remap[name]]['xsec']
                scale = tgt_lumi * xsec / events
                obj.Scale(scale)
            if 'data_obs' in name:
                tgt_lumi = sample_cfg[remap['data_obs_%s' % chn]]['lumi']
                obj.SetTitle('lumi:%.1f fb^{-1} (13 TeV, %s)' % (tgt_lumi / 1000., year))
        for sample in samples:
            if 'data_obs' not in sample:
                tgt_lumi = 35916.764758
                events = sample_cfg[remap[sample]]['events']
                xsec = sample_cfg[remap[sample]]['xsec']
                scale = tgt_lumi * xsec / events
                print '%-50s %-20.1f %-12.2f %-12.2f' % (remap[sample], events, xsec, scale)

    for chn in ['e', 'm']:
        for S in ['barrel_%s' % chn, 'endcap_%s' % chn]:
            ratio = hists[chn]['%s_iso_l_sig_t' % S]['p0_pt']['W_J'].Clone()
            ratio.Divide(hists[chn]['%s_iso_l_sig_l' % S]['p0_pt']['W_J'])
            hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_RATIO'] = ratio.Clone()
            ratio.Multiply(hists[chn]['%s_iso_t_sig_l' % S]['p0_pt']['W_J'])
            hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_EST'] = ratio
            hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_SDB'] = hists[chn]['%s_iso_t_sig_l' % S]['p0_pt']['W_J'].Clone()
            true_ratio = hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_J'].Clone()
            true_ratio.Divide(hists[chn]['%s_iso_t_sig_l' % S]['p0_pt']['W_J'])
            hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_TRUERATIO'] = true_ratio.Clone()

            bias_ratio = hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_RATIO'].Clone()
            bias_ratio.Divide(hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_TRUERATIO'])
            hists[chn]['%s_iso_t_sig_t' % S]['p0_pt']['W_BIASRATIO'] = bias_ratio.Clone()

else:
    fin = ROOT.TFile('fakestudy_2016_photon_fakes.root')
    hists = Node()
    TDirToNode(fin, node=hists)
    fin.Close()

def Analyse2D(hist, xvar, yvar, xcut, ycut, xstart=None, ystart=None):
    xcutbin = hist.GetXaxis().FindFixBin(xcut)
    ycutbin = hist.GetYaxis().FindFixBin(ycut)
    xbins = hist.GetNbinsX()
    ybins = hist.GetNbinsY()
    xstartbin = hist.GetXaxis().FindFixBin(xcut if xstart is None else xstart)
    ystartbin = hist.GetYaxis().FindFixBin(ycut if ystart is None else ystart)
    if xstart is None:
        xstartbin += 1
    if ystart is None:
        ystartbin += 1

    print 'X: found bin [%f,%f] for cut value %f of variable %s' % (hist.GetXaxis().GetBinLowEdge(xcutbin), hist.GetXaxis().GetBinUpEdge(xcutbin), xcut, xvar)
    print 'Y: found bin [%f,%f] for cut value %f of variable %s' % (hist.GetYaxis().GetBinLowEdge(ycutbin), hist.GetYaxis().GetBinUpEdge(ycutbin), ycut, yvar)

    print 'X: SR = [%f, %f], CR = [%f, ...]' % (hist.GetXaxis().GetBinLowEdge(1), hist.GetXaxis().GetBinUpEdge(xcutbin), hist.GetXaxis().GetBinLowEdge(xstartbin))
    print 'Y: SR = [%f, %f], CR = [%f, ...]' % (hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(ycutbin), hist.GetYaxis().GetBinLowEdge(ystartbin))
    h_bias = hist.Clone()
    h_bias.Reset()

    h_stat = hist.Clone()
    h_stat.Reset()

    y_A = hist.Integral(1, xcutbin, 1, ycutbin)
    for ix in xrange(xstartbin, xbins + 1):
        for iy in xrange(ystartbin, ybins + 1):
            y_B = hist.Integral(xstartbin, ix, 1, ycutbin)
            y_C = hist.Integral(1, xcutbin, ystartbin, iy)
            y_D = hist.Integral(xstartbin, ix, ystartbin, iy)
            if y_A > 0 and y_B > 0 and y_C > 0 and y_D > 0:
                h_bias.SetBinContent(ix, iy, (y_C / y_D) / (y_A / y_B))
            # h_stat.SetBinContent(ix, iy, (y_C / y_D) / (y_A / y_B))
    return (h_bias, h_stat)
fout = ROOT.TFile('fakestudy_%s_%s.root' % (year, args.task), 'RECREATE')

for ivar in xrange(0, len(p0_vars)):
    for jvar in xrange(ivar + 1, len(p0_vars)):
        xname = p0_vars[ivar]
        yname = p0_vars[jvar]
        label = '%s_vs_%s' % (xname, yname)
        if xname == 'p0_chiso' and yname == 'p0_sigma':
            h_bias, h_stat = Analyse2D(hists['2D']['m']['full_barrel_m'][label]['W'], xname, yname, p0_cuts[xname], p0_cuts[yname], xstart=4.0, ystart=0.01100)
        else:
            h_bias, h_stat = Analyse2D(hists['2D']['m']['full_barrel_m'][label]['W'], xname, yname, p0_cuts[xname], p0_cuts[yname])
        hists['2D']['m']['full_barrel_m'][label]['BIAS'] = h_bias
        hists['2D']['m']['full_barrel_m'][label]['STAT'] = h_stat

for path, node in hists.ListNodes(withObjects=True):
    print path
    if path.startswith('2D'):
        continue

    # node['VV'] = HistSum(dibosons)
    # node['TTG'] = HistSum(ttgammas)
    # if year in ['2017', '2018']:
    #     node['TT'] = HistSum(tt_samples)

    # for P in ['F', 'E', 'J', 'R']:
    #     node['VV_%s' % P] = HistSum(['%s_%s' % (lbl, P) for lbl in dibosons])
    #     node['TTG_%s' % P] = HistSum(['%s_%s' % (lbl, P) for lbl in ttgammas])
    #     if year in ['2017', '2018']:
    #         node['TT_%s' % P] = HistSum(['%s_%s' % (lbl, P) for lbl in tt_samples])

    # node['Total_R'] = HistSum(['WG_NLO', 'TT_R', 'DY_XZG_R', 'ZG_IZG_R', 'VV_R'])
    # node['Total_E'] = HistSum(['W_E', 'TT_E', 'DY_E', 'VV_E'])
    # node['Total_J'] = HistSum(['W_J', 'TT_J', 'DY_J', 'VV_J'])
    # node['Total_F'] = HistSum(['W_F', 'TT_F', 'DY_F', 'VV_F'])

NodeToTDir(fout, hists)

fout.Close()
