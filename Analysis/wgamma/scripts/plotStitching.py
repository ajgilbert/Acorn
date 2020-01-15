import ROOT
import argparse
import json
import itertools
from copy import deepcopy
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
from Acorn.Analysis.plottemplates import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
# plot.ModTDRStyle(width=900, height=900)

samples = {
    'LowZG': 'ZGToLLG-lowMLL-amcatnloFXFX',
    'ZG': 'ZGToLLG-amcatnloFXFX-stitched'
}
hists = Node()

sample_files = {
    key: 'root://eoscms.cern.ch//store/user/agilbert/191216-full/wgamma_2016_v4/WGamma_' + value + '.root' for (key, value) in samples.iteritems()
}

for sa in samples:
    hists['ptl'][sa] = Hist('TH1D', sample=sa, var=['ptl'], binning=(20, 0, 100), sel='1', wt='wt')
    hists['mll_sf'][sa] = Hist('TH1D', sample=sa, var=['mll_sf'], binning=(20, 0, 200), sel='1', wt='wt')

# hists = dict()
# f_inc = ROOT.TFile('nlo_lhe_pt_inc.root')
# hists['inc'] = f_inc.Get('lhe_pt')
# f_120 = ROOT.TFile('nlo_lhe_pt_130_after_PS.root')
# hists['120'] = f_120.Get('lhe_pt')

# tot = hists['inc'].Integral(0, hists['inc'].GetNbinsX() + 1)
# hi = hists['inc'].Integral(hists['inc'].GetXaxis().FindFixBin(120), hists['inc'].GetNbinsX() + 1)

# print hi / tot

# hists['inc'].Scale(1. / tot)
# hists['120'].Scale((hi / tot) / hists['120'].Integral(0, hists['120'].GetNbinsX() + 1))
MultiDraw(hists, sample_files, 'SampleStitching', mt_cores=4, mt_thresh=1E5)

plotcfg = DEFAULT_CFG
plotcfg.update({
    'type': 'multihist',
    'legend_pos': [0.55, 0.75, 0.90, 0.91],
    'ratio_pad_frac': 0.33,
    'top_title_right': '',
    'logy': True,
    'logy_min': 1E-6
    })

x_titles = {
    'ptl': ['LHE lepton p_{T}', 'GeV'],
    'mll_sf': ['LHE m_{ll}', 'GeV']
}

for var in ['ptl', 'mll_sf']:
    hist_dict = {}
    for opath, objname, obj in hists[var].ListObjects(depth=0):
        hist_dict[objname] = obj

    MakeMultiHistPlot('zg_stitching_%s' % var,
                        outdir='.',
                        hists=hist_dict,
                        cfg=UpdateDict(plotcfg, {
                            'x_title': x_titles[var],
                            'ratio_y_title': 'Ratio',
                            'top_title_right': '',
                            'ratio_y_range': [0.5, 1.5],
                            }),
                        layout=[
                            {'name': 'LowZG', 'legend': 'Low m_{ll} only'},
                            {'name': 'ZG', 'legend': 'Merged sample'},
                        ],
                        ratios=[
                            {'num': 'LowZG', 'den': 'ZG'},
                        ])



f_baseline = ROOT.TFile('output_2016_baseline_191210-WWG-split.root')

hist_dict = {}
for sa in ['WWTo1L1Nu2Q_R', 'WWG_R', 'WWTo1L1Nu2Q_IWWG_R']:
    hist_dict[sa] = f_baseline.Get('/m/baseline_m_mZ_veto/p0_pt/%s' % sa)


MakeMultiHistPlot('wwg_stitching',
                    outdir='.',
                    hists=hist_dict,
                    cfg=UpdateDict(plotcfg, {
                        'x_title': ['Photon p_{T}', 'GeV'],
                        'ratio_y_title': 'Ratio',
                        'ratio': False,
                        'logy': False,
                        'top_title_right': '',
                        'ratio_y_range': [0.5, 1.5],
                        'rebinvar': [30, 50, 80, 100, 200]
                        }),
                    layout=[
                        {'name': 'WWG_R', 'legend': 'WW#gamma'},
                        {'name': 'WWTo1L1Nu2Q_R', 'legend': 'WW (W#rightarrowqq FSR only)'},
                        {'name': 'WWTo1L1Nu2Q_IWWG_R', 'legend': 'WW (remaining #gamma FSR/ISR)'},
                    ])



f_baseline = ROOT.TFile('output_2016_baseline_191216.root')

hist_dict = {}
for sa in ['TGJets_R', 'ST_t_antitop_R', 'ST_t_top_R']:
    hist_dict[sa] = f_baseline.Get('/m/baseline_m_mZ_veto/p0_pt/%s' % sa)


MakeMultiHistPlot('tg_stitching',
                    outdir='.',
                    hists=hist_dict,
                    cfg=UpdateDict(plotcfg, {
                        'x_title': ['Photon p_{T}', 'GeV'],
                        'ratio_y_title': 'Ratio',
                        'ratio': False,
                        'logy': False,
                        'top_title_right': '',
                        'ratio_y_range': [0.5, 1.5],
                        'rebinvar': [30, 50, 80, 100, 200]
                        }),
                    layout=[
                        {'name': 'TGJets_R', 'legend': 't#gamma (t-channel)'},
                        {'name': 'ST_R', 'entries': ['ST_t_antitop_R', 'ST_t_top_R'], 'legend': 't+jets (t-channel)'},
                    ])

# # Draw efficiencies for the stitched sample
# for year, var in itertools.product(years, ['l0_pt', 'p0_pt']):
#     hist_dict = {}
#     for opath, objname, obj in hists[year][var]['WG'].ListObjects(depth=0):
#         hist_dict[objname] = obj

#     # Draw N/N-1 efficiencies
#     MakeMultiHistPlot('efficiencies_m_%s_%s' % (year, var),
#                       outdir=outdir,
#                       hists=hist_dict,
#                       cfg=UpdateDict(plotcfg, {
#                           'x_title': x_titles[var],
#                           'ratio_y_title': '#varepsilon(N/N-1)',
#                           'top_title_right': '%.1f fb^{-1} (13 TeV, %s)' % (lumi_fb[year], year),
#                           'ratio_y_range': [0.55, 1.05],
#                           }),
#                       layout=[
#                           {'name': 'baseline_m', 'legend': 'Baseline'},
#                           {'name': 'muon_id', 'legend': ' + Muon ID'},
#                           {'name': 'muon_iso', 'legend': ' + Muon Iso'},
#                           {'name': 'muon_trg', 'legend': ' + Muon Trg(24/27) OR (50)'},
#                           {'name': 'muon_trg_2', 'legend': ' + Muon Trg(24/27)'},
#                           {'name': 'photon_id_m', 'legend': ' + Photon ID'},
#                           {'name': 'photon_pix_m', 'legend': ' + Photon pix veto'},
#                       ],
#                       ratios=[
#                           {'num': 'muon_id', 'den': 'baseline_m', 'type': 'binomial'},
#                           {'num': 'muon_iso', 'den': 'muon_id', 'type': 'binomial'},
#                           {'num': 'muon_trg', 'den': 'muon_iso', 'type': 'binomial'},
#                           {'num': 'muon_trg_2', 'den': 'muon_trg', 'type': 'binomial'},
#                           {'num': 'photon_id_m', 'den': 'muon_trg_2', 'type': 'binomial'},
#                           {'num': 'photon_pix_m', 'den': 'photon_id_m', 'type': 'binomial'},
#                       ]
#     )

