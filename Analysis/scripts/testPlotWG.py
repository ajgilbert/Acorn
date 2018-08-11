import ROOT
import json
from pprint import pprint
from collections import defaultdict
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)


tname = 'WGDataAnalysis'

samples = {
    'DY': 'output/090818/wgamma_2016_v2/WGamma/DYJetsToLL_M-50-madgraphMLM.root',
    'data_obs': 'output/090818/wgamma_2016_v2/WGamma/SingleMuon.root',
    'TT': 'output/090818/wgamma_2016_v2/WGamma/TT-powheg.root',
    'WG': 'output/090818/wgamma_2016_v2/WGamma/WGToLNuG-madgraphMLM.root',
    'W': 'output/090818/wgamma_2016_v2/WGamma/WJetsToLNu-madgraphMLM.root'
}

remap = {
    'DY': 'DYJetsToLL_M-50-madgraphMLM',
    'data_obs': 'SingleMuon',
    'TT': 'TT-powheg',
    'WG': 'WGToLNuG-madgraphMLM',
    'W': 'WJetsToLNu-madgraphMLM'
}

hists = Node()

# for sample in samples:
#     hists['mt'][sample] = Hist('TH1F', (20, 0., 150.), sample, ['mt'], sel='trg_m')
    # hists['eta_1'][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_1'], sel='pt_1 > %s && trg_1' % ptcut)
    # hists['eta_2'][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_2'], sel='pt_1 > %s && trg_1' % ptcut)
    # hists['pt_1'][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_1'], sel='pt_1 > %s && trg_1' % ptcut)
    # hists['pt_2'][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_2'], sel='pt_1 > %s && trg_1' % ptcut)
for sample in samples:
    hists['mt_weighted'][sample] = Hist('TH1D', (20, 0., 250.), sample, ['m0met_mt'], sel='m0_trg && n_m==1', wt='wt_pu*wt_m0*wt_trg_m0')
    # hists['pt_met_weighted'][sample] = Hist('TH1F', (20, 0., 150.), sample, ['pt_met'], sel='trg_m', wt='wt_pu*wt_m*wt_trg_m')
    # hists['pt_m_weighted'][sample] = Hist('TH1F', [0,20,40,60,80,100,125,150,175,200,250,300,400,500], sample, ['pt_m'], sel='trg_m', wt='wt_pu*wt_m*wt_trg_m')
    # hists['pt_p_weighted'][sample] = Hist('TH1F', [0,20,40,60,80,100,125,150,175,200,250,300,400,500], sample, ['pt_p'], sel='trg_m', wt='wt_pu*wt_m*wt_trg_m')
    # hists['eta_1_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_1'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')
    # hists['eta_2_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_2'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')
    # hists['pt_1_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_1'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')
    # hists['pt_2_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_2'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')


MultiDraw(hists, samples, tname)

with open('input/cfg_wgamma_2016_v2.json') as jsonfile:
    cfg = json.load(jsonfile)
    sample_cfg = cfg['samples']

for sample in samples:
    f = ROOT.TFile(samples[sample])
    sample_cfg[remap[sample]]['events'] = f.Get('counters').GetBinContent(2)
    f.Close()

for path, name, obj in hists.ListObjects():
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
