import CombineHarvester.CombineTools.plotting as plot
import ROOT
import argparse
import sys
import os
import fnmatch
from copy import deepcopy

ROOT.TH1.SetDefaultSumw2(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()


def GetListOfDirectories(file, dirlist=None, curr=[]):
    if dirlist is None:
        dirlist = list()
        file.cd()
    ROOT.gDirectory.cd('/' + '/'.join(curr))
    dirs_to_proc = []
    has_other_objects = False
    for key in ROOT.gDirectory.GetListOfKeys():
        if key.GetClassName() == 'TDirectoryFile':
            dirs_to_proc.append(key.GetName())
        else:
            has_other_objects = True
    if has_other_objects:
        # This directory should be added to the list
        dirlist.append('/' + '/'.join(curr))
    for dirname in dirs_to_proc:
        GetListOfDirectories(file, dirlist, curr + [dirname])
    return dirlist


def GetHistsInDir(file, dir):
    res = {}
    file.cd()
    ROOT.gDirectory.cd(dir)
    for key in ROOT.gDirectory.GetListOfKeys():
        if key.GetClassName() != 'TDirectoryFile':
            obj = key.ReadObj()
            res[obj.GetName()] = obj
    return res


LAYOUTS = {
    "wgamma": [
        ('VV', {
            'entries': ['VVTo2L2Nu', 'WWTo1L1Nu2Q', 'WZTo1L1Nu2Q', 'WZTo1L3Nu', 'WZTo2L2Q', 'WZTo3LNu', 'ZZTo2L2Q', 'ZZTo4L'],
            'legend': 'Diboson',
            'color': ROOT.TColor.GetColor(248, 206, 104)
        }
        ),
        ('DY', {
            'entries': ['DY_R'],
            'legend': 'Z#rightarrowll',
            'color': ROOT.TColor.GetColor(100, 192, 232)
        }
        ),
        ('DY_F', {
            'entries': ['DY_F'],
            'legend': 'Z#rightarrowll (F)',
            'color': 2
        }
        ),
        ('TT', {
            'entries': ['TT_R'],
            'legend': 't#bar{t}',
            'color': ROOT.TColor.GetColor(155, 152, 204)
        }
        ),
        ('TT_F', {
            'entries': ['TT_F'],
            'legend': 't#bar{t}',
            'color': 4
        }
        ),
        ('W', {
            'entries': ['W_F'],
            'legend': 'W#rightarrowl#nu',
            'color': ROOT.TColor.GetColor(222, 90, 106)
        }
        ),
        ('WG', {
            # 'entries': ['WG'],
            'entries': ['WG_p_acc', 'WG_n_acc'],
            'legend': 'W#rightarrowl#nu+#gamma (in)',
            'color': ROOT.TColor.GetColor(119, 213, 217)

        }
        ),
        ('WG_OOA', {
            'entries': ['WG_p_ooa', 'WG_n_ooa'],
            'legend': 'W#rightarrowl#nu+#gamma (out)',
            'color': 12
        }
        ),
        # ('ZTT', {
        #     'entries': ['ZTT'],
        #     'legend': 'Z#rightarrow#tau#tau',
        #     'color': ROOT.TColor.GetColor(248, 206, 104)
        # }
        # )
    ]
}


config_by_setting = {
    "xtitle": {
        '*/m0met_mt': ('m_{T}(#mu,p_{T}^{miss})', 'GeV'),
        '*/m0_pt': ('Muon p_{T}', 'GeV'),
        '*/m0_eta': ('Muon #eta', ''),
        '*/m0_iso': ('Muon Iso', 'GeV'),
        '*/m1_pt': ('Subleading muon p_{T}', 'GeV'),
        '*/m1_eta': ('Subleading muon #eta', ''),
        '*/m0m1_M': ('m_{#mu^{+}#mu^{-}}', 'GeV'),
        '*/m0m1_dr': ('#DeltaR(#mu^{+},#mu^{-})', ''),
        '*/met': ('p_{T}^{miss}', 'GeV'),
        '*/p0_pt': ('Photon p_{T}', 'GeV'),
        '*/p0_eta': ('Photon #eta', ''),
        '*/m0p0_dr': ('#DeltaR(#mu,#gamma)', ''),
        '*/m0p0_M': ('m_{#mu#gamma}', 'GeV'),
        '*/p0_chiso': ('Photon I_{charged}', 'GeV'),
        '*/p0_neiso': ('Photon I_{neutral}', 'GeV'),
        '*/p0_phiso': ('Photon I_{photon}', 'GeV'),
        '*/p0_hovere': ('Photon H/E', ''),
        '*/p0_sigma': ('Photon #sigma_{I#etaI#eta}^{full 5x5}', ''),
        '*/p0_haspix': ('Photon hasPixelSeed', ''),
        '*/abs(reco_phi)': ('Reconstructed #phi', ''),
    }
}

variants_by_path = {
    "*": {
    },
    "*/p0_pt": {
        "prefix": "logy_",
        "logy": True
    },
    "*/p0_pt": {
        "prefix": "zoom_",
        "range": (0, 200)
    }
}

# Derive fake template as A = B * (C / D) from data


def Subtracted(hists, target, subtract, zero=True):
    res = hists[target].Clone()
    for sub in subtract:
        res.Add(hists[sub], -1.)
    for i in xrange(1, res.GetNbinsX() + 1):
        if res.GetBinContent(i) < 0.:
            res.SetBinContent(i, 0.)
    return res


def SafeDivide(h1, h2):
    res = h1.Clone()
    for i in xrange(1, res.GetNbinsX() + 1):
        if h2.GetBinContent(i) == 0.:
            res.SetBinContent(i, 0.)
        else:
            res.SetBinContent(i, h1.GetBinContent(i) / h2.GetBinContent(i))
    return res


def Multiply(h1, h2):
    res = h1.Clone()
    for i in xrange(1, res.GetNbinsX() + 1):
        res.SetBinContent(i, h1.GetBinContent(i) * h2.GetBinContent(i))
    return res


def DoPhotonFakes(file, b, c, d):
    h_b = GetHistsInDir(file, b)
    h_b['data_sub'] = Subtracted(h_b, 'data_obs', ['DY', 'WG'])
    h_c = GetHistsInDir(file, c)
    h_c['data_sub'] = Subtracted(h_c, 'data_obs', ['DY', 'WG'])
    h_d = GetHistsInDir(file, d)
    h_d['data_sub'] = Subtracted(h_d, 'data_obs', ['DY', 'WG'])
    h_div = SafeDivide(h_c['data_sub'], h_d['data_sub'])
    h_res = Multiply(h_b['data_sub'], h_div)
    return h_res


def MakePlot(name, outdir, hists, cfg):
    copyhists = {}
    for hname, h in hists.iteritems():
        copyhists[hname] = h.Clone()
        if cfg['divwidth']:
            copyhists[hname].Scale(1., 'width')

    hists = copyhists
    # Canvas and pads
    canv = ROOT.TCanvas(name, name)
    pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

    # Get the data and create axis hist
    h_data = hists['data_obs']
    # h_data = Getter(file, '%s/data_obs' % target, True)
    if isinstance(h_data, ROOT.TH2):
        print 'TH2: aborting!'
        return

    h_axes = [h_data.Clone() for x in pads]
    for h in h_axes:
        if len(cfg['range']):
            h.GetXaxis().SetRangeUser(*cfg['range'])
        h.Reset()

    build_h_tot = True
    h_tot = None
    if 'TotalProcs' in hists:
        h_tot = hists['TotalProcs']
        build_h_tot = False

    x_title = cfg['xtitle'][0]
    units = cfg['xtitle'][1]

    if x_title == '' and h_data.GetXaxis().GetTitle() != '':
        x_title = h_data.GetXaxis().GetTitle()

    if ':' in x_title:
        units = x_title.split(':')[1]
        x_title = x_title.split(':')[0]

    if cfg['logy']:
        pads[0].SetLogy()
        h_axes[0].SetMinimum(0.1)
    plot.StandardAxes(h_axes[1].GetXaxis(), h_axes[0].GetYaxis(), x_title, units)
    h_axes[0].Draw()

    # A dict to keep track of the hists
    h_store = {}

    layout = LAYOUTS['wgamma']

    stack = ROOT.THStack()
    legend = ROOT.TLegend(0.67, 0.86 - 0.04*len(layout), 0.90, 0.91, '', 'NBNDC')

    # h_tot = None

    for ele in layout:
        info = ele[1]
        hist = hists[info['entries'][0]]
        plot.Set(hist, FillColor=info['color'], Title=info['legend'])
        if len(info['entries']) > 1:
            for other in info['entries'][1:]:
                hist.Add(hists[other])
        h_store[ele[0]] = hist
        if build_h_tot:
            if h_tot is None:
                h_tot = hist.Clone()
            else:
                h_tot.Add(hist)
        stack.Add(hist)

    h_tot.SetFillColor(plot.CreateTransparentColor(12, 0.3))
    h_tot.SetMarkerSize(0)

    legend.AddEntry(h_data, 'Observed', 'PL')
    for ele in reversed(layout):
        legend.AddEntry(h_store[ele[0]], '', 'F')
    bkg_uncert_label = 'Stat. Uncertainty'
    if not build_h_tot:
        bkg_uncert_label = 'Uncertainty'
    legend.AddEntry(h_tot, bkg_uncert_label, 'F')

    stack.Draw('HISTSAME')
    h_tot.Draw("E2SAME")
    h_data.Draw('E0SAME')

    # if args.logy:
    #     h_axes[0].SetMinimum(args.y_min)
    #     pads[0].SetLogy(True)

    plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.35)
    legend.Draw()
    # plot.FixBoxPadding(pads[0], legend, 0.05)

    # Do the ratio plot
    pads[1].cd()
    pads[1].SetGrid(0, 1)
    h_axes[1].Draw()

    r_data = plot.MakeRatioHist(h_data, h_tot, True, False)
    r_tot = plot.MakeRatioHist(h_tot, h_tot, True, False)
    r_tot.Draw('E2SAME')
    r_data.Draw('SAME')

    plot.SetupTwoPadSplitAsRatio(
        pads, plot.GetAxisHist(
            pads[0]), plot.GetAxisHist(pads[1]), 'Obs/Exp', True, 0.61, 1.69)

    # Go back and tidy up the axes and frame
    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    # CMS logo
    plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)
    plot.DrawTitle(pads[0], '35.9 fb^{-1} (13 TeV)', 3)

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

args = parser.parse_args()

default_cfg = {
    'prefix': '',
    'postfix': '',
    'logy': False,
    'logx': False,
    'range': [],
    'rebin': 0,
    'rebinvar': [],
    'xtitle': '',
    'ytitle': 'Events',
    'divwidth': True
}

filename, dirfilter = args.input.split(':')
print filename
file = ROOT.TFile(filename)

target_list = GetListOfDirectories(file)
print target_list

filtered_list = [x for x in target_list if fnmatch.fnmatch(x, dirfilter)]
print filtered_list

draw_list = []

for path in filtered_list:
    # for now work on the assumption that the last component of the path will be the actual filename
    split_path = path.split('/')[:-1]
    name = path.split('/')[-1]
    target_dir = os.path.join(args.output, *split_path)
    os.system('mkdir -p %s' % target_dir)
    hists = GetHistsInDir(file, path)
    # if 'w_highmt_pho' in path:
    #     hists['W_DAT'] = DoPhotonFakes(file, 'w_hmt_pho_iso_l_sig_t/%s' % name, 'w_hmt_pho_iso_t_sig_l/%s' % name, 'w_hmt_pho_iso_l_sig_l/%s' % name)

    # print hists

    plotcfg = dict(default_cfg)
    for setting, vardict in config_by_setting.iteritems():
        for pathkey, val in vardict.iteritems():
            if fnmatch.fnmatch(path, pathkey):
                plotcfg[setting] = val
                print 'Path %s, setting %s, to value %s' % (path, setting, val)
    for pathkey, varcfg in variants_by_path.iteritems():
        if fnmatch.fnmatch(path, pathkey):
            varplotcfg = dict(plotcfg)
            varplotcfg.update(varcfg)
            MakePlot(name, target_dir, hists, varplotcfg)


