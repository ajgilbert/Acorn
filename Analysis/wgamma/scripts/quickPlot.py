import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
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


def HistSum(hdict, label_list):
    return sum([hdict[X] for X in label_list[1:]], hdict[label_list[0]])


default_cfg = {
    # full filename will be [outdir]/[prefix]name[postfix].[ext]:
    'outdir': '',               # output directory
    'prefix': '',               # filename prefix
    'postfix': '',              # filename postfix
    'logx': False,              # Draw x-axis in log-scale
    'logy': False,              # Draw y-axis in log-scale
    'ratio': True,             # Draw the ratio plot?
    'fraction': False,             # Draw the ratio plot?
    'ratio_y_range': [0.61, 1.39],  # Range of the ratio y-axis
    'x_range': [],                  # Restrict the x-axis range shown
    'rebin': 0,                     # Rebin by this factor
    'rebinvar': [],                 # Rebin to this list of bin edges
    'x_title': '',               # Title on the x-axis
    'y_title': 'Events',         # Title on the y-axis
    'divwidth': True,           # Divide all histogram contents by their bin widths
    'layout': 'data_fakes',
    'legend_pos': [0.67, 0.66, 0.90, 0.91], # Legend position in NDC (x1, y1, x2, y2)
    'legend_padding': 0.05,      # Automatically increase the y-axis range to ensure the legend does not overlap the histograms (argument is fraction of frame height to pad)
    'data_name': 'data_obs',     # Name of the TH1 to take for data
    'main_logo': 'CMS',
    'sub_logo': 'Internal',
    'top_title_right': '137.0 fb^{-1} (13 TeV)',
    'top_title_left': '',
    'hide_data': False,
    'auto_top_title_right': True,
    'overlays': [
        # {
        #     "name": "systUp",
        #     "entries": "total",
        #     "hist_postfix": "_no_p0_wt",
        #     "legend": "No Photon ID SF",
        #     "color": 2
        # },
        # {
        #     "name": "systDown",
        #     "entries": "total",
        #     "hist_postfix": "_CMS_ele_fake_pDown",
        #     "legend": "systDown",
        #     "color": 4
        # }
    ]
}

config_by_setting = {
    "x_title": {
        '*/n_vtx': ('Number of vertices', ''),
        '*/n_veto_m': ('Number of veto muons', ''),
        '*/n_veto_e': ('Number of veto electrons', ''),
        'm/*/l0met_mt': ('m_{T}(#mu,p_{T}^{miss})', 'GeV'),
        'e/*/l0met_mt': ('m_{T}(e,p_{T}^{miss})', 'GeV'),
        'm/*/l0_pt': ('Muon p_{T}', 'GeV'),
        'e/*/l0_pt': ('Electron p_{T}', 'GeV'),
        'm/*/l0_eta': ('Muon #eta', ''),
        'e/*/l0_eta': ('Electron #eta', ''),
        'm/*/l0_phi': ('Muon #phi', ''),
        'e/*/l0_phi': ('Electron #phi', ''),
        'm/*/l0_iso': ('Muon Iso', 'GeV'),
        'e/*/l0_iso': ('Electron Iso', 'GeV'),
        'm/*/l1_pt': ('Subleading muon p_{T}', 'GeV'),
        'm/*/l1_eta': ('Subleading muon #eta', ''),
        'm/*/l0l1_M': ('m_{#mu^{+}#mu^{-}}', 'GeV'),
        'm/*/l0l1_pt': ('p_{T}^{#mu^{+}#mu^{-}}', 'GeV'),
        'm/*/l0l1_dr': ('#DeltaR(#mu^{+},#mu^{-})', ''),
        'e/*/l1_pt': ('Subleading electron p_{T}', 'GeV'),
        'e/*/l1_eta': ('Subleading electron #eta', ''),
        'e/*/l0l1_M': ('m_{e^{+}e^{-}}', 'GeV'),
        'e/*/l0l1_pt': ('p_{T}^{e^{+}e^{-}}', 'GeV'),
        'e/*/l0l1_dr': ('#DeltaR(e^{+},e^{-})', ''),
        '*/*/l0j0_dphi': ('#Delta#phi(l,jet)', ''),
        '*/*/l0met_dphi': ('#Delta#phi(l,MET)', ''),
        '*/*/j0_pt': ('Leading jet p_{T}', 'GeV'),
        '*/met': ('p_{T}^{miss}', 'GeV'),
        '*/met_phi': ('p_{T}^{miss} #phi', ''),
        '*/tk_met': ('Track p_{T}^{miss}', 'GeV'),
        '*/tk_met_phi': ('Track p_{T}^{miss} #phi', ''),
        '*/puppi_met': ('PUPPI p_{T}^{miss}', 'GeV'),
        '*/puppi_met_phi': ('PUPPI p_{T}^{miss} #phi', ''),
        '*/gen_p0_pt': ('Gen. Photon p_{T}', 'GeV'),
        '*/p0_pt': ('Photon p_{T}', 'GeV'),
        '*/p0_eta': ('Photon #eta', ''),
        '*/p0_phi': ('Photon #phi', ''),
        '*/p0_truth': ('Photon truth match', ''),
        'm/*/l0p0_dr': ('#DeltaR(#mu,#gamma)', ''),
        'm/*/l0p0_M': ('m_{#mu#gamma}', 'GeV'),
        'e/*/l0p0_dr': ('#DeltaR(e,#gamma)', ''),
        'e/*/l0p0_M': ('m_{e#gamma}', 'GeV'),
        '*/p0_chiso': ('Photon I_{charged}', 'GeV'),
        '*/p0_neiso': ('Photon I_{neutral}', 'GeV'),
        '*/p0_phiso': ('Photon I_{photon}', 'GeV'),
        '*/p0_hovere': ('Photon H/E', ''),
        '*/p0_sigma': ('Photon #sigma_{I#etaI#eta}^{full 5x5}', ''),
        '*/p0_haspix': ('Photon hasPixelSeed', ''),
        '*/p0_eveto': ('Photon passElectronVeto', ''),
        '*/abs(reco_phi)': ('Reconstructed #phi', ''),
        '*/abs(gen_phi)': ('Gen. #phi', ''),
        '*/abs(true_phi)': ('True #phi', ''),
        '*/abs(reco_phi_f)': ('Reconstructed #phi_{f}', ''),
        '*/abs(reco_puppi_phi_f)': ('Reconstructed #phi_{f}', ''),
        '*/abs(gen_phi_f)': ('Gen. #phi_{f}', ''),
        '*/abs(true_phi_f)': ('True #phi_{f}', ''),
        '*/wt_def': ('Default weight', ''),
        '*/wt_pu': ('Pileup weight', ''),
        '*/wt_l0': ('Lepton weight', ''),
        '*/wt_l1': ('Sub-leading lepton weight', ''),
        '*/wt_trg_l0': ('Lepton trigger weight', ''),
        '*/wt_p0': ('Photon weight', ''),
        '*/wt_p0_e_fake': ('e#rightarrow#gamma weight', ''),
        '*/1': ('1', ''),
        '*/0.5': ('0.5', '')
    },
    "layout": {
        'm/cr_Zmm/*': 'pure_mc_zll',
        'e/cr_Zee/*': 'pure_mc_zll',
        'm/lepton_fakes_m/*': 'pure_mc',
        'e/lepton_fakes_e/*': 'pure_mc',
    }
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
    ("*/p0_chiso", {
            "prefix": "logy_",
            "logy": True}),
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
    # ("*/p0_pt", {
    #         "prefix": "fr_barrel_",
    #         "rebinvar": [30, 35, 40, 50, 60, 80, 100, 300],
    #         "logy": True})
]


def MakePlot(name, outdir, hists, cfg, layouts):
    copyhists = {}
    for hname, h in hists.iteritems():
        if len(cfg['rebinvar']):
            copyhists[hname] = VariableRebin(h, cfg['rebinvar'])
        elif cfg['rebin'] > 0:
            copyhists[hname] = h.Clone()
            copyhists[hname].Rebin(cfg['rebin'])
        else:
            copyhists[hname] = h.Clone()
        if cfg['divwidth']:
            copyhists[hname].Scale(1., 'width')

    hists = copyhists

    # Canvas and pads
    canv = ROOT.TCanvas(name, name)
    if cfg['ratio'] or cfg['fraction']:
        pads = plot.TwoPadSplit(0.27, 0.01, 0.01)
    else:
        pads = plot.OnePad()

    # Get the data and create axis hist
    h_data = hists[cfg['data_name']]
    # h_data = Getter(file, '%s/data_obs' % target, True)
    if isinstance(h_data, ROOT.TH2):
        print 'TH2: aborting!'
        return

    h_axes = [h_data.Clone() for x in pads]
    for h in h_axes:
        if len(cfg['x_range']):
            h.GetXaxis().SetRangeUser(*cfg['x_range'])
        h.Reset()

    build_h_tot = True
    h_tot = None
    if 'TotalProcs' in hists:
        h_tot = hists['TotalProcs']
        build_h_tot = False

    x_title = cfg['x_title'][0]
    units = cfg['x_title'][1]

    if x_title == '' and h_data.GetXaxis().GetTitle() != '':
        x_title = h_data.GetXaxis().GetTitle()

    if ':' in x_title:
        units = x_title.split(':')[1]
        x_title = x_title.split(':')[0]

    if cfg['logy']:
        pads[0].SetLogy()
        h_axes[0].SetMinimum(0.001)

    if cfg['ratio'] or cfg['fraction']:
        plot.StandardAxes(h_axes[1].GetXaxis(), h_axes[0].GetYaxis(), x_title, units)
    else:
        plot.StandardAxes(h_axes[0].GetXaxis(), h_axes[0].GetYaxis(), x_title, units)
    h_axes[0].Draw()

    # A dict to keep track of the hists
    h_store = {}

    layout = layouts[cfg['layout']]

    stack = ROOT.THStack()
    legend = ROOT.TLegend(*(cfg['legend_pos'] + ['', 'NBNDC']))

    all_input_hists = []
    for info in layout:
        all_input_hists.extend(info['entries'])
        hist = hists[info['entries'][0]].Clone()
        col = info['color']
        if isinstance(col, list):
            col = ROOT.TColor.GetColor(*col)
        plot.Set(hist, FillColor=col, Title=info['legend'])
        if len(info['entries']) > 1:
            for other in info['entries'][1:]:
                hist.Add(hists[other])
        h_store[info['name']] = hist
        if build_h_tot:
            if h_tot is None:
                h_tot = hist.Clone()
            else:
                h_tot.Add(hist)
        stack.Add(hist)

    h_tot.SetFillColor(plot.CreateTransparentColor(12, 0.3))
    h_tot.SetMarkerSize(0)

    # Build overlays
    for info in cfg['overlays']:
        hist = None
        input_list = []
        if isinstance(info['entries'], str):
            input_list = list(all_input_hists)
        else:
            input_list = list(info['entries'])
        updated_list = []
        for xh in input_list:
            if xh + info['hist_postfix'] in hists:
                updated_list.append(xh + info['hist_postfix'])
            else:
                updated_list.append(xh)
        print updated_list
        hist = HistSum(hists, updated_list)
        col = info['color']
        if isinstance(col, list):
            col = ROOT.TColor.GetColor(*col)
        plot.Set(hist, LineColor=col, LineWidth=2, MarkerSize=0, Title=info['legend'])
        for ib in xrange(1, hist.GetNbinsX() + 1):
            hist.SetBinError(ib, 1E-7)
        h_store[info['name']] = hist


    legend.AddEntry(h_data, 'Observed', 'PL')
    for ele in reversed(layout):
        legend.AddEntry(h_store[ele['name']], '', 'F')
    bkg_uncert_label = 'Stat. Uncertainty'
    if not build_h_tot:
        bkg_uncert_label = 'Uncertainty'
    legend.AddEntry(h_tot, bkg_uncert_label, 'F')

    for info in cfg['overlays']:
        legend.AddEntry(h_store[info['name']], info['legend'], 'L')

    stack.Draw('HISTSAME')
    h_tot.Draw("E2SAME")

    for info in cfg['overlays']:
        h_store[info['name']].Draw('SAME')

    if not cfg['hide_data']:
        h_data.Draw('E0SAME')

    plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.35)
    legend.Draw()
    if cfg['legend_padding'] > 0.:
        print h_axes[0].GetMinimum(), h_axes[0].GetMaximum()
        if not h_axes[0].GetMaximum() == 0.:
            plot.FixBoxPadding(pads[0], legend, cfg['legend_padding'])

    # Do the ratio plot
    r_store = {}
    if cfg['ratio'] or cfg['fraction']:
        pads[1].cd()
        pads[1].SetGrid(0, 1)
        h_axes[1].Draw()

        if cfg['ratio']:
            r_data = plot.MakeRatioHist(h_data, h_tot, True, False)
            r_tot = plot.MakeRatioHist(h_tot, h_tot, True, False)
            r_tot.Draw('E2SAME')
            if not cfg['hide_data']:
                r_data.Draw('SAME')
            for info in cfg['overlays']:
                r_store[info['name']] = plot.MakeRatioHist(h_store[info['name']], h_tot, False, False)
                r_store[info['name']].Draw('SAME')

            plot.SetupTwoPadSplitAsRatio(
                pads, plot.GetAxisHist(
                    pads[0]), plot.GetAxisHist(pads[1]), 'Obs/Exp', True, cfg['ratio_y_range'][0], cfg['ratio_y_range'][1])
        if cfg['fraction']:
            r_frac = plot.MakeRatioHist(h_tot, h_data, True, True)
            r_frac.Draw('SAME')
            plot.SetupTwoPadSplitAsRatio(
                pads, plot.GetAxisHist(
                    pads[0]), plot.GetAxisHist(pads[1]), 'Exp/Obs', True, 0.0, 0.5)

    # Go back and tidy up the axes and frame
    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    # CMS logo
    plot.DrawCMSLogo(pads[0], cfg['main_logo'], cfg['sub_logo'], 11, 0.045, 0.05, 1.0, '', 1.0)
    plot.DrawTitle(pads[0], cfg['top_title_left'], 1)
    if cfg['auto_top_title_right']:
        title_right = h_data.GetTitle()
        if title_right.startswith('lumi:'):
            plot.DrawTitle(pads[0], title_right.replace('lumi:', ''), 3)
    else:
        plot.DrawTitle(pads[0], cfg['top_title_right'], 3)

    latex = ROOT.TLatex()
    plot.Set(latex, NDC=None, TextFont=42, TextSize=0.03)
    latex.DrawLatex(0.20, 0.75, args.title)
    # plot.DrawTitle(pads[0], args.title, 1)

    # ... and we're done
    canv.Print(outdir + '/' + cfg['prefix'] + name + cfg['postfix'] + '.png')
    canv.Print(outdir + '/' + cfg['prefix'] + name + cfg['postfix'] + '.pdf')

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

default_cfg['layout'] = args.layout

if args.lumi is not None:
    default_cfg['top_title_right'] = args.lumi
    default_cfg['auto_top_title_right'] = False

with open(args.layout_file) as jsonfile:
    layouts = json.load(jsonfile)

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
        for pathkey, val in vardict.iteritems():
            if fnmatch.fnmatch(path, pathkey):
                plotcfg[setting] = val
                print 'Path %s, setting %s, to value %s' % (path, setting, val)
    for pathkey, varcfg in variants_by_path:
        if fnmatch.fnmatch(path, pathkey):
            varplotcfg = dict(plotcfg)
            varplotcfg.update(varcfg)
            MakePlot(name, target_dir, hists, varplotcfg, layouts)
