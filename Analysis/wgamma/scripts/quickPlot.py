import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from Acorn.Analysis.plottemplates import *
import ROOT
import argparse
import json
import os
import fnmatch
from copy import deepcopy
from array import array

ROOT.TH1.SetDefaultSumw2(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()


# for ix in xrange(50):
#     plot.CreateTransparentColor(12, 0.0 + 0.001 * ix)
# def HistSum(hdict, label_list):
#     return sum([hdict[X] for X in label_list[1:]], hdict[label_list[0]])


mod_cfg = {
    'hide_data': False,
    'global_hist_opts': {
        'draw_opts': 'HIST',
        'legend_opts': 'F',
        'marker_size': 0.6,
        'line_width': 1
    }
}

config_by_setting = {
    "x_title": [
        ('*/n_vtx', ('Number of vertices', '')),
        ('*/n_veto_m', ('Number of veto muons', '')),
        ('*/n_veto_e', ('Number of veto electrons', '')),
        ('*/n_alt_veto_m', ('Number of veto muons', '')),
        ('*/n_alt_veto_e', ('Number of veto electrons', '')),
        ('m/*/l0met_mt', ('m_{T}(#mu,p_{T}^{miss})', 'GeV')),
        ('e/*/l0met_mt', ('m_{T}(e,p_{T}^{miss})', 'GeV')),
        ('m/*/l0_pt', ('Muon p_{T}', 'GeV')),
        ('e/*/l0_pt', ('Electron p_{T}', 'GeV')),
        ('m/*/l0_eta', ('Muon #eta', '')),
        ('e/*/l0_eta', ('Electron #eta', '')),
        ('m/*/l0_phi', ('Muon #phi', '')),
        ('e/*/l0_phi', ('Electron #phi', '')),
        ('m/*/l0_iso', ('Muon Iso', 'GeV')),
        ('e/*/l0_iso', ('Electron Iso', 'GeV')),
        ('m/*/l1_pt', ('Subleading muon p_{T}', 'GeV')),
        ('m/*/l1_eta', ('Subleading muon #eta', '')),
        ('m/*/l0l1_M', ('m_{#mu^{+}#mu^{-}}', 'GeV')),
        ('m/*/l0l1_pt', ('p_{T}^{#mu^{+}#mu^{-}}', 'GeV')),
        ('m/*/l0l1_dr', ('#DeltaR(#mu^{+},#mu^{-})', '')),
        ('e/*/l1_pt', ('Subleading electron p_{T}', 'GeV')),
        ('e/*/l1_eta', ('Subleading electron #eta', '')),
        ('e/*/l0l1_M', ('m_{e^{+}e^{-}}', 'GeV')),
        ('e/*/l0l1_pt', ('p_{T}^{e^{+}e^{-}}', 'GeV')),
        ('e/*/l0l1_dr', ('#DeltaR(e^{+},e^{-})', '')),
        ('*/*/l0j0_dphi', ('#Delta#phi(l,jet)', '')),
        ('*/*/l0met_dphi', ('#Delta#phi(l,MET)', '')),
        ('*/*/j0_pt', ('Leading jet p_{T}', 'GeV')),
        ('*/met', ('p_{T}^{miss}', 'GeV')),
        ('*/met_phi', ('p_{T}^{miss} #phi', '')),
        ('*/tk_met', ('Track p_{T}^{miss}', 'GeV')),
        ('*/tk_met_phi', ('Track p_{T}^{miss} #phi', '')),
        ('*/puppi_met', ('PUPPI p_{T}^{miss}', 'GeV')),
        ('*/puppi_met_phi', ('PUPPI p_{T}^{miss} #phi', '')),
        ('*/gen_p0_pt', ('Gen. Photon p_{T}', 'GeV')),
        ('*/p0_pt', ('Photon p_{T}', 'GeV')),
        ('*/p0_eta', ('Photon #eta', '')),
        ('*/p0_phi', ('Photon #phi', '')),
        ('*/p0_truth', ('Photon truth match', '')),
        ('m/*/l0p0_dr', ('#DeltaR(#mu,#gamma)', '')),
        ('*/*/l0p0_deta', ('#Delta#eta(l,#gamma)', '')),
        ('m/*/l0p0_M', ('m_{#mu#gamma}', 'GeV')),
        ('e/*/l0p0_dr', ('#DeltaR(e,#gamma)', '')),
        ('e/*/l0p0_M', ('m_{e#gamma}', 'GeV')),
        ('*/*/mt_cluster*', ('m_{T}^{cluster}', 'GeV')),
        ('*/p0_worstiso', ('Photon I_{worst}', 'GeV')),
        ('*/p0_chiso', ('Photon I_{charged}', 'GeV')),
        ('*/p0_neiso', ('Photon I_{neutral}', 'GeV')),
        ('*/p0_phiso', ('Photon I_{photon}', 'GeV')),
        ('*/p0_hovere', ('Photon H/E', '')),
        ('*/p0_sigma', ('Photon #sigma_{I#etaI#eta}^{full 5x5}', '')),
        ('*/p0_haspix', ('Photon hasPixelSeed', '')),
        ('*/p0_eveto', ('Photon passElectronVeto', '')),
        ('*/abs(reco_phi)', ('Reconstructed #phi', '')),
        ('*/abs(gen_phi)', ('Gen. #phi', '')),
        ('*/abs(true_phi)', ('True #phi', '')),
        ('*/abs(reco_phi_f)', ('Reconstructed #phi_{f}', '')),
        ('*/abs(reco_puppi_phi_f)', ('Reconstructed #phi_{f}', '')),
        ('*/abs(gen_phi_f)', ('Gen. #phi_{f}', '')),
        ('*/abs(true_phi_f)', ('True #phi_{f}', '')),
        ('*/wt_def', ('Default weight', '')),
        ('*/wt_pu', ('Pileup weight', '')),
        ('*/wt_l0', ('Lepton weight', '')),
        ('*/wt_l1', ('Sub-leading lepton weight', '')),
        ('*/wt_trg_l0', ('Lepton trigger weight', '')),
        ('*/wt_p0', ('Photon weight', '')),
        ('*/wt_p0_e_fake', ('e#rightarrow#gamma weight', '')),
        ('*/gen_mll', ('Gen. m_{ll}', 'GeV')),
        ('*/n_all_j', ('Number of jets', '')),
        ('*/n_all_btag_j', ('Number of b-tagged jets', '')),
        ('*/1', ('1', '')),
        ('*/0.5', ('0.5', ''))
    ],
    "data_name": [
        ('*/wt_*', 'WG_R')
    ],
    "hide_data": [
        ('*/wt_*', True)
    ],
    "layout": [
        ('m/cr_Zmm/*', 'pure_mc_zll'),
        ('e/cr_Zee/*', 'pure_mc_zll'),
        ('m/lepton_fakes_m/*', 'pure_mc'),
        ('e/lepton_fakes_e/*', 'pure_mc'),
        ('e/*/wt_*', 'data_fakes_e_NODATA'),
        ('m/*/wt_*', 'data_fakes_m_NODATA'),
    ]
}

variants_by_path = [
    ("*", {}),
    # ("*", {
    #         "postfix": "_mc",
    #         "layout": "pure_mc"
    #     }),
    ("*/p0_pt", {
            "prefix": "logy_",
            "logy": True}),
    ("*/gen_mll", {
            "prefix": "mask_",
            "hide_data": True
        }),
    ("*/p0_chiso", {
            "prefix": "logy_",
            "logy": True}),
    ("*/p0_worstiso", {
            "prefix": "logy_",
            "logy": True,
            "rebin": 2}),
    ("*/p0_pt", {
            "prefix": "zoom_",
            "x_range": (0, 200)}),
    ("*/l0l1_M", {
            "prefix": "rebinned_",
            "rebin": 2,
            "logy": True}),
    ("*/p0_pt", {
            "prefix": "zoom_logy_",
            "x_range": (0, 200),
            "logy": True}),
    ("*/p0_eta", {
            "prefix": "wide_",
            "rebinvar": [-3.0, -2.7, -1.5, -0.9, 0, 0.9, 1.5, 2.7, 3.0]}),
    ("*/l0p0_deta", {
            "prefix": "fid_",
            "rebinvar": [-5, -3.4, -2.6, -1.8, -1.4, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.4, 1.8, 2.6, 3.4, 5]}),
]


parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='Output of PostFitShapes or PostFitShapesFromWorkspace, specified as FILE:BIN')
parser.add_argument('--output', '-o', default='.', help='top level output directory')
parser.add_argument('--channel', '-c', default='mt', choices=['mt', 'et', 'em', 'tt', 'mm', 'generic', 'wgamma'], help='Channel')
parser.add_argument('--x-title', default='', help='x-axis variable, without GeV')
parser.add_argument('--logy', action='store_true')
parser.add_argument('--y-min', type=float, default=1)
parser.add_argument('--title', default='')
parser.add_argument('--lumi', default=None)
parser.add_argument('--layout-file', '-l', default='layouts.json')
parser.add_argument('--layout', default='data_fakes')

args = parser.parse_args()

default_cfg = dict(DEFAULT_CFG)
default_cfg.update(mod_cfg)
default_cfg['layout'] = args.layout

if args.lumi is not None:
    default_cfg['top_title_right'] = args.lumi
    default_cfg['auto_top_title_right'] = False

with open(args.layout_file) as jsonfile:
    layouts = json.load(jsonfile)
    for l in layouts.keys():
        layouts[l + "_NODATA"] = [X for X in layouts[l] if (not 'data' in X['name'])]
filename, dirfilter = args.input.split(':')
print filename
file = ROOT.TFile(filename)

node = TDirToNode(file)

made_dirs = set()

for path, subnode in node.ListNodes(withObjects=True):
    print path
    if not fnmatch.fnmatch(path, dirfilter):
        continue
    # for now work on the assumption that the last component of the path will be the actual filename
    split_path = path.split('/')[:-1]
    name = path.split('/')[-1]
    target_dir = os.path.join(args.output, *split_path)
    if target_dir not in made_dirs:
        os.system('mkdir -p %s' % target_dir)
        made_dirs.add(target_dir)
    hists = {}
    for opath, objname, obj in subnode.ListObjects(depth=0):
        hists[objname] = obj

    plotcfg = dict(default_cfg)
    for setting, vardict in config_by_setting.iteritems():
        for pathkey, val in vardict:
            if fnmatch.fnmatch(path, pathkey):
                plotcfg[setting] = val
                print 'Path %s, setting %s, to value %s' % (path, setting, val)
    for pathkey, varcfg in variants_by_path:
        if fnmatch.fnmatch(path, pathkey):
            varplotcfg = dict(plotcfg)
            varplotcfg.update(varcfg)
            MakeMultiHistPlot(name, target_dir, hists, varplotcfg, layouts[varplotcfg['layout']])
            # MakePlot(name, target_dir, hists, varplotcfg, layouts)
