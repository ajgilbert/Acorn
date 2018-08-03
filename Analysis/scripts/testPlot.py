import ROOT
import json
from pprint import pprint
from collections import defaultdict
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)


tname = 'DiMuonAnalysis'

samples = {
    'DY': 'output/060618/wgamma_2016_v1/Main/DYJetsToLL_M-50-madgraphMLM.root',
    'data_obs': 'output/060618/wgamma_2016_v1/Main/SingleMuon.root',
    'TT': 'output/060618/wgamma_2016_v1/Main/TT-powheg.root',
    'WG': 'output/060618/wgamma_2016_v1/Main/WGToLNuG-madgraphMLM.root'
}

remap = {
    'DY': 'DYJetsToLL_M-50-madgraphMLM',
    'data_obs': 'SingleMuon',
    'TT': 'TT-powheg',
    'WG': 'WGToLNuG-madgraphMLM'
}

hists = Node()

for ptcut in range(30, 100, 10):
    for sample in samples:
        hists['m_ll']['pt_%s' % ptcut][sample] = Hist('TH1F', (60, 60., 120.), sample, ['m_ll'], sel='pt_1 > %s && trg_1' % ptcut)
        hists['eta_1']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_1'], sel='pt_1 > %s && trg_1' % ptcut)
        hists['eta_2']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_2'], sel='pt_1 > %s && trg_1' % ptcut)
        hists['pt_1']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_1'], sel='pt_1 > %s && trg_1' % ptcut)
        hists['pt_2']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_2'], sel='pt_1 > %s && trg_1' % ptcut)
    for sample in samples:
        hists['m_ll_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (60, 60., 120.), sample, ['m_ll'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')
        hists['eta_1_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_1'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')
        hists['eta_2_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, -2.4, 2.4), sample, ['eta_2'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')
        hists['pt_1_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_1'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')
        hists['pt_2_weighted']['pt_%s' % ptcut][sample] = Hist('TH1F', (40, 0, 100), sample, ['pt_2'], sel='pt_1 > %s && trg_1' % ptcut, wt='wt_pu*wt_1*wt_2*wt_trg1')


MultiDraw(hists, samples, tname)

with open('input/cfg_wgamma_2016_v1.json') as jsonfile:
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
