import ROOT
import argparse
from copy import deepcopy
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
from Acorn.Analysis.plottemplates import *

# import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
# ROOT.TH1.SetDefaultSumw2(True)
plot.ModTDRStyle()


parser = argparse.ArgumentParser()
# parser.add_argument('--input', '-i', help='Input directory')
# parser.add_argument('--output', '-o', default='.', help='top level output directory')
# parser.add_argument('--channel', '-c', default='mt', choices=['mt', 'et', 'em', 'tt', 'mm', 'generic', 'wgamma'], help='Channel')
# parser.add_argument('--x-title', default='', help='x-axis variable, without GeV')
# parser.add_argument('--logy', action='store_true')
# parser.add_argument('--y-min', type=float, default=1)
# parser.add_argument('--title', default='')
# parser.add_argument('--layout-file', '-l', default='layouts.json')
# parser.add_argument('--layout', default='data_fakes')

args = parser.parse_args()

samples = {
    'WG': 'root://eoscms.cern.ch//store/cmst3/user/agilbert/190205/wgamma_2017_v3/WGamma/WGToLNuG-madgraphMLM-stitched.root'
    # 'WG': 'output/test/wgamma_2017_v3/WGamma/WGToLNuG-madgraphMLM-PtG-500_0.root'
    # 'WG': 'output/test/wgamma_2017_v3/WGamma/WGToLNuG-madgraphMLM_0.root'
}

# remap = {
#     'WG': 'WGToLNuG-EFT_pTA_300_inf-madgraphMLM'
# }

# fout = ROOT.TFile('%s_w_%s.root' % (args.output, pm_label), 'RECREATE')
hists = Node()

list_of_vars = [
    ('lhe_p0_pt', 'lhe_p0_pt', (50, 0, 1000)),
    ('p0_pt', 'gen_p0_pt', [0, 10, 20, 30, 40, 50, 100, 150, 210, 300, 420, 600, 850, 1200]),
    ('p0_eta', 'gen_p0_eta', (20, -3, 3)),
    ('l0_pt', 'gen_l0_pt', [0, 10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500]),
    ('l0_eta', 'gen_l0_eta', (20, -3, 3)),
]

for label, var, binning in list_of_vars:
    X = SelectionManager()
    X['baseline'] = 'gen_pdgid==13 && is_wg_gen && gen_l0_pt > 30 && abs(gen_l0_eta) < 2.4 && gen_p0_pt > 30 && abs(gen_p0_eta) < 2.5 && (abs(gen_p0_eta) < 1.4442 || abs(gen_p0_eta) > 1.566) && gen_l0p0_dr > 0.7'
    X['muon_id'] = '$baseline && n_pre_m>=1 && l0_pdgid == 13'
    X['muon_iso'] = '$muon_id && l0_iso < 0.15'
    X['muon_trg'] = '$muon_iso && (l0_trg || (l0_pt > 55 && l0_trg_2))'
    X['muon_trg_2'] = '$muon_iso && l0_trg'
    X['photon_pre'] = '$muon_trg_2 && n_pre_p>=1'
    X['photon_id'] = '$photon_pre && p0_medium'
    X['photon_pix'] = '$photon_id && !p0_haspix'

    for sel in ['baseline', 'muon_id', 'muon_iso', 'muon_trg', 'muon_trg_2', 'photon_pre', 'photon_id', 'photon_pix']:
        hists[label][sel] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.get('$' + sel), wt='wt_pu*wt_def')

# X.get('$photon_pix', printlevel=1)

MultiDraw(hists, samples, 'WGDataAnalysis')

plotcfg = DEFAULT_CFG
plotcfg.update({
    'type': 'multihist',
    'ratio_y_range': [0.65, 1.05],
    'legend_pos': [0.55, 0.66, 0.90, 0.91],
    'ratio_pad_frac': 0.33,
    'ratio_y_title': '#varepsilon(N/N-1)',
    'top_title_right': '41.5 fb^{-1} (13 TeV, 2017)',
    'logy': True
    })

x_titles = {
    'l0_pt': ['Gen. lepton p_{T}', 'GeV'],
    'p0_pt': ['Gen. photon p_{T}', 'GeV'],
    'lhe_p0_pt': ['LHE photon p_{T}', 'GeV'],
}

for var in ['l0_pt', 'p0_pt']:
    hist_dict = {}
    for opath, objname, obj in hists[var].ListObjects(depth=0):
        hist_dict[objname] = obj

    cfg = deepcopy(plotcfg)
    cfg.update({
        'x_title': x_titles[var]
        })

    hist_list = [
            {'name': 'baseline', 'legend': 'Baseline', 'color': 2},
            {'name': 'muon_id', 'legend': ' + Muon ID', 'color': 1},
            {'name': 'muon_iso', 'legend': ' + Muon Iso', 'color': 8},
            {'name': 'muon_trg', 'legend': ' + Muon Trg(24/27) OR (50)', 'color': 5},
            {'name': 'muon_trg_2', 'legend': ' + Muon Trg(24/27)', 'color': 4}
        ]
    MakeMultiHistPlot('efficiencies_2017_%s' % var, '.', hist_dict, cfg, hist_list, ratios=[
            {'num': 'muon_id', 'den': 'baseline', 'type': 'binomial'},
            {'num': 'muon_iso', 'den': 'muon_id', 'type': 'binomial'},
            {'num': 'muon_trg', 'den': 'muon_iso', 'type': 'binomial'},
            {'num': 'muon_trg_2', 'den': 'muon_trg', 'type': 'binomial'},
            # {'num': 'photon_ge1', 'den': 'muon_trg'},
            # {'num': 'photon_pre', 'den': 'muon_id'},
            # {'num': 'photon_id', 'den': 'photon_pre'},
            # {'num': 'photon_pix', 'den': 'photon_id'},
        ])
    cfg.update({
        'ratio_y_title': '#varepsilon(N/baseline)',
        'ratio_y_range': [0.35, 1.05],
        })
    MakeMultiHistPlot('efficiencies_cumulative_2017_%s' % var, '.', hist_dict, cfg, hist_list, ratios=[
            {'num': 'muon_id', 'den': 'baseline', 'type': 'binomial'},
            {'num': 'muon_iso', 'den': 'baseline', 'type': 'binomial'},
            {'num': 'muon_trg', 'den': 'baseline', 'type': 'binomial'},
            {'num': 'muon_trg_2', 'den': 'baseline', 'type': 'binomial'},
            # {'num': 'photon_ge1', 'den': 'muon_trg'},
            # {'num': 'photon_pre', 'den': 'muon_id'},
            # {'num': 'photon_id', 'den': 'photon_pre'},
            # {'num': 'photon_pix', 'den': 'photon_id'},
        ])

for var in ['lhe_p0_pt']:
    hist_dict = {}
    for opath, objname, obj in hists[var].ListObjects(depth=0):
        hist_dict[objname] = obj

    cfg = deepcopy(plotcfg)
    cfg.update({
        'x_title': x_titles[var]
        })

    hist_list = [
            {'name': 'baseline', 'legend': 'Baseline', 'color': 2}
        ]
    MakeMultiHistPlot('lhe_2017_%s' % var, '.', hist_dict, cfg, hist_list)
