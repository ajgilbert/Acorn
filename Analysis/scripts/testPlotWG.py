import ROOT
import json
import sys
from pprint import pprint
from collections import defaultdict
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)


tname = 'WGDataAnalysis'

prefix = 'output/130818-reduced/wgamma_2016_v2/WGamma/'

samples = {
    'DY': 'DYJetsToLL_M-50-madgraphMLM.root',
    'data_obs': 'SingleMuon.root',
    'TT': 'TT-powheg.root',
    'WG': 'WGToLNuG-madgraphMLM.root',
    'WG-St': 'WGToLNuG-madgraphMLM-stitched.root',
    'WG-PtG-130': 'WGToLNuG-madgraphMLM-PtG-130.root',
    'WG-PtG-500': 'WGToLNuG-madgraphMLM-PtG-500.root',
    'W': 'WJetsToLNu-madgraphMLM.root',
    'VVTo2L2Nu': 'VVTo2L2Nu-amcatnloFXFX.root',
    'WWTo1L1Nu2Q': 'WWTo1L1Nu2Q-amcatnloFXFX.root',
    'WZTo1L1Nu2Q': 'WZTo1L1Nu2Q-amcatnloFXFX.root',
    'WZTo1L3Nu': 'WZTo1L3Nu-amcatnloFXFX.root',
    'WZTo2L2Q': 'WZTo2L2Q-amcatnloFXFX.root',
    'WZTo3LNu': 'WZTo3LNu-amcatnloFXFX.root',
    'ZZTo2L2Q': 'ZZTo2L2Q-amcatnloFXFX.root',
    'ZZTo4L': 'ZZTo4L-amcatnloFXFX.root',
}

for sa in samples:
    samples[sa] = (prefix + samples[sa])

remap = {
    'DY': 'DYJetsToLL_M-50-madgraphMLM',
    'data_obs': 'SingleMuon',
    'TT': 'TT-powheg',
    'WG': 'WGToLNuG-madgraphMLM',
    'WG-St': 'WGToLNuG-madgraphMLM-stitched',
    'WG-PtG-130': 'WGToLNuG-madgraphMLM-PtG-130',
    'WG-PtG-500': 'WGToLNuG-madgraphMLM-PtG-500',
    'W': 'WJetsToLNu-madgraphMLM',
    'VVTo2L2Nu': 'VVTo2L2Nu-amcatnloFXFX',
    'WWTo1L1Nu2Q': 'WWTo1L1Nu2Q-amcatnloFXFX',
    'WZTo1L1Nu2Q': 'WZTo1L1Nu2Q-amcatnloFXFX',
    'WZTo1L3Nu': 'WZTo1L3Nu-amcatnloFXFX',
    'WZTo2L2Q': 'WZTo2L2Q-amcatnloFXFX',
    'WZTo3LNu': 'WZTo3LNu-amcatnloFXFX',
    'ZZTo2L2Q': 'ZZTo2L2Q-amcatnloFXFX',
    'ZZTo4L': 'ZZTo4L-amcatnloFXFX'
}

hists = Node()

X = SelectionManager()
X.Set('baseline', sel='m0_trg && n_m >= 1', wt='wt_def*wt_pu*wt_m0*wt_trg_m0')
X.Derive('w_inc', base='baseline', sel='n_m==1')
X.Derive('w_inc_iso_m', base='w_inc', sel='m0_iso < 0.15')
# X.Derive('w_inc_aiso_m', base='w_inc', sel='m0_iso > 0.2 && m0_iso < 0.5')
# X.Derive('z_inc', 'baseline', sel='n_m==2 && m0m1_os && m0m1_dr>0.5', wt='wt_m1')

# X.Derive('w_lowmt', 'w_inc_iso_m', sel='m0met_mt<30')
# X.Derive('w_lowmt_aiso_m', 'w_inc_aiso_m', sel='m0met_mt<30')

# X.Derive('w_lowmt_pho_rec', 'w_lowmt', sel='n_p==1 && m0p0_dr>0.7')
# X.Derive('w_lowmt_pho_rec_aiso_m', 'w_lowmt_aiso_m', sel='n_p==1 && m0p0_dr>0.7')

X.Derive('w_highmt', 'w_inc_iso_m', sel='met>40')
# X.Derive('w_highmt_aiso_m', 'w_inc_aiso_m', sel='m0met_mt>60 && met>40')

X.Derive('w_highmt_pho_rec', 'w_highmt', sel='n_p==1 && m0p0_dr>0.7')
# X.Derive('w_highmt_pho_rec', 'w_highmt', sel='n_p==1 && m0p0_dr>0.7 && abs(p0_eta) < 1.4442')
# X.Derive('w_highmt_pho_rec_aiso_m', 'w_highmt_aiso_m', sel='n_p==1 && m0p0_dr>0.7')

# X.Derive('w_highmt_pho', 'w_highmt_pho_rec', sel='p0_medium && !p0_haspix', wt='wt_p0')
X.Derive('w_highmt_pho', 'w_highmt_pho_rec', sel='p0_medium && !p0_haspix && m0_pt>80 && met>80 && p0_pt>150 && m0p0_dr>3.0', wt='wt_p0')
X.Derive('w_highmt_pho_b1', 'w_highmt_pho', sel='p0_pt>150 && p0_pt<210', wt='wt_p0')
X.Derive('w_highmt_pho_b2', 'w_highmt_pho', sel='p0_pt>210 && p0_pt<300', wt='wt_p0')
X.Derive('w_highmt_pho_b3', 'w_highmt_pho', sel='p0_pt>300 && p0_pt<420', wt='wt_p0')
X.Derive('w_highmt_pho_b4', 'w_highmt_pho', sel='p0_pt>420 && p0_pt<1200', wt='wt_p0')
# X.Derive('w_highmt_pho_b4', 'w_highmt_pho', sel='p0_pt>420 && p0_pt<600', wt='wt_p0')
# X.Derive('w_highmt_pho_b4', 'w_highmt_pho', sel='p0_pt>600 && p0_pt<850', wt='wt_p0')
# X.Derive('w_highmt_pho_b5', 'w_highmt_pho', sel='p0_pt>850 && p0_pt<1200', wt='wt_p0')
# X.Derive('w_highmt_pho', 'w_highmt_pho_rec', sel='p0_medium && !p0_haspix && m0_pt>80 && met>80 && p0_pt>300 && m0p0_dr>3.0', wt='wt_p0')

# X.Derive('w_highmt_pho_aiso_m', 'w_highmt_pho_rec_aiso_m', sel='p0_medium && !p0_haspix', wt='wt_p0')


X.Derive('w_hmt_pho_pre', 'w_highmt_pho_rec', sel='p0_medium_noch && !p0_haspix', wt='wt_p0')
X.Derive('w_hmt_pho_iso_l', 'w_hmt_pho_pre', sel='p0_chiso > 2.5 && p0_chiso < 7', wt='wt_p0')
X.Derive('w_hmt_pho_sig_l', 'w_hmt_pho_pre', sel='p0_sigma > 0.01022', wt='wt_p0')
X.Derive('w_hmt_pho_iso_l_sig_t', 'w_hmt_pho_pre', sel='p0_chiso > 2.5 && p0_chiso < 7 && p0_sigma < 0.01022', wt='wt_p0')
X.Derive('w_hmt_pho_iso_l_sig_l', 'w_hmt_pho_pre', sel='p0_chiso > 2.5 && p0_chiso < 7 && p0_sigma > 0.01022', wt='wt_p0')
X.Derive('w_hmt_pho_iso_t_sig_l', 'w_hmt_pho_pre', sel='p0_chiso < 0.441 && p0_sigma > 0.01022', wt='wt_p0')


drawvars = [
    # ('m0met_mt', (30, 0., 200.)),
    ('m0_pt', (40, 0., 150.)),
    # ('m0_eta', (20, -3.0, 3.0)),
    # ('m0_iso', (40, 0, 2.0)),
    # ('m1_pt', (40, 0., 150.)),
    # ('m1_eta', (20, -3.0, 3.0)),
    # ('m0m1_M', (40, 60, 120)),
    # ('m0m1_dr', (20, 0., 5.)),
    ('met', (20, 0., 200.)),
    ('p0_pt', [0, 10, 20, 30, 40, 50, 60, 80, 100, 150, 210, 300, 420, 600, 850, 1200]),
    # ('p0_eta', (20, -3.0, 3.0)),
    ('m0p0_dr', (20, 0., 5.)),
    # ('m0p0_M', (40, 60, 120)),
    # ('p0_chiso', (40, 0, 20.0)),
    # ('p0_neiso', (40, 0, 20.0)),
    # ('p0_phiso', (40, 0, 20.0)),
    # ('p0_hovere', (20, 0., 0.5)),
    # ('p0_sigma', (40, 0., 0.050)),
    # ('p0_haspix', (2, -0.5, 1.5)),
    ('abs(reco_phi)', (5, 0, 3.15)),
]

for sel in X.storage.keys():
    for var, binning in drawvars:
        for sample in samples:
            hists[sel][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel=X.sel('$'+sel), wt=X.wt('$'+sel))
        hists[sel][var]['W_R'] = Hist('TH1D', sample='W', var=[var], binning=binning, sel=X.sel('$'+sel + ' && p0_isprompt'), wt=X.wt('$'+sel))
        hists[sel][var]['W_F'] = Hist('TH1D', sample='W', var=[var], binning=binning, sel=X.sel('$'+sel + ' && !p0_isprompt'), wt=X.wt('$'+sel))
        hists[sel][var]['DY_R'] = Hist('TH1D', sample='DY', var=[var], binning=binning, sel=X.sel('$'+sel + ' && p0_isprompt'), wt=X.wt('$'+sel))
        hists[sel][var]['DY_F'] = Hist('TH1D', sample='DY', var=[var], binning=binning, sel=X.sel('$'+sel + ' && !p0_isprompt'), wt=X.wt('$'+sel))
        hists[sel][var]['TT_R'] = Hist('TH1D', sample='TT', var=[var], binning=binning, sel=X.sel('$'+sel + ' && p0_isprompt'), wt=X.wt('$'+sel))
        hists[sel][var]['TT_F'] = Hist('TH1D', sample='TT', var=[var], binning=binning, sel=X.sel('$'+sel + ' && !p0_isprompt'), wt=X.wt('$'+sel))

MultiDraw(hists, samples, tname)

with open('input/cfg_wgamma_2016_v2.json') as jsonfile:
    cfg = json.load(jsonfile)
    sample_cfg = cfg['samples']

for sample in samples:
    f = ROOT.TFile(samples[sample])
    sample_cfg[remap[sample]]['events'] = f.Get('counters').GetBinContent(2)
    f.Close()

for path, hname, obj in hists.ListObjects():
    name = obj.sample
    if name is not 'data_obs':
        tgt_lumi = sample_cfg[remap['data_obs']]['lumi']
        events = sample_cfg[remap[name]]['events']
        xsec = sample_cfg[remap[name]]['xsec']
        scale = tgt_lumi * xsec / events
        obj.Scale(scale)

fout = ROOT.TFile('output_2016.root', 'RECREATE')

for path, name, obj in hists.ListObjects():
    WriteToTFile(obj, fout, path, name)

fout.Close()
