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
X.Derive('z_inc', 'baseline', sel='n_m==2 && m0m1_os && m0m1_dr>0.5', wt='wt_m1')
X.Derive('w_highmt', 'w_inc', sel='m0met_mt>60 && met>40')
X.Derive('w_highmt_pho_rec', 'w_highmt', sel='n_p==1 && m0p0_dr>0.7')
X.Derive('w_highmt_pho', 'w_highmt_pho_rec', sel='p0_medium && !p0_haspix', wt='wt_p0')
X.Derive('w_highmt_invch', 'w_highmt_pho_rec', sel='p0_medium_noch && p0_chiso > 2. && p0_chiso < 8.')


drawvars = [
    ('m0met_mt', (30, 0., 200.)),
    ('m0_pt', (40, 0., 150.)),
    ('m0_eta', (20, -3.0, 3.0)),
    ('m1_pt', (40, 0., 150.)),
    ('m1_eta', (20, -3.0, 3.0)),
    ('m0m1_M', (40, 60, 120)),
    ('m0m1_dr', (20, 0., 5.)),
    ('met', (20, 0., 200.)),
    ('p0_pt', [0, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 600, 1000, 2000]),
    ('p0_eta', (20, -3.0, 3.0)),
    ('m0p0_dr', (20, 0., 5.)),
    ('m0p0_M', (40, 60, 120)),
    ('p0_chiso', (40, 0, 20.0)),
    ('p0_neiso', (40, 0, 20.0)),
    ('p0_phiso', (40, 0, 20.0)),
    ('p0_hovere', (20, 0., 0.5)),
    ('p0_sigma', (40, 0., 0.050)),
    ('p0_haspix', (2, -0.5, 1.5)),
]

for sel in ['baseline', 'w_inc', 'z_inc', 'w_highmt', 'w_highmt_pho_rec', 'w_highmt_pho', 'w_highmt_invch']:
    for var, binning in drawvars:
        for sample in samples:
            hists[sel][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel=X.sel('$'+sel), wt=X.wt('$'+sel))
        hists[sel][var]['W_R'] = Hist('TH1D', sample='W', var=[var], binning=binning, sel=X.sel('$'+sel + ' && p0_isprompt'), wt=X.wt('$'+sel))
        hists[sel][var]['W_F'] = Hist('TH1D', sample='W', var=[var], binning=binning, sel=X.sel('$'+sel + ' && !p0_isprompt'), wt=X.wt('$'+sel))

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
