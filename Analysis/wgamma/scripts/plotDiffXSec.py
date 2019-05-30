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

parser = argparse.ArgumentParser()
parser.add_argument('--selection', default='eft_region')
parser.add_argument('--scheme', default='phi_f_binned')
parser.add_argument('--years', default='2016,2017,2018')
parser.add_argument('--ratio', action='store_true')
parser.add_argument('--charge', default='p')
args = parser.parse_args()


plot.ModTDRStyle(width=1000)
ROOT.gStyle.SetEndErrorSize(0)

def DoScaleEnvelope(node, nominal):
    h = node[nominal].Clone()
    h_alts = []
    for i in range(6):
        h_alts.append(node[nominal + '__sc_%i' % i])
    for ix in xrange(1, h.GetNbinsX() + 1):
        for iy in xrange(1, h.GetNbinsY() + 1):
            max_dev = max([abs(x.GetBinContent(ix, iy) - h.GetBinContent(ix, iy)) for x in h_alts])
            h.SetBinError(ix, iy, max_dev)
    return h

def DivideGraphByHist(gr, hist):
    res = gr.Clone()
    for i in xrange(gr.GetN()):
        res.GetY()[i] = res.GetY()[i] / hist.GetBinContent(i + 1)
    if type(res) is ROOT.TGraphAsymmErrors:
        for i in xrange(gr.GetN()):
            res.GetEYhigh()[i] = res.GetEYhigh()[i]/hist.GetBinContent(i + 1)
            res.GetEYlow()[i] = res.GetEYlow()[i]/hist.GetBinContent(i + 1)
    return res


def SetupPads(split_points, gaps_left, gaps_right, ratio=None):
    pads = []
    ratio_pads = []
    l_margin = ROOT.gStyle.GetPadLeftMargin()
    r_margin = ROOT.gStyle.GetPadRightMargin()
    t_margin = ROOT.gStyle.GetPadTopMargin()
    b_margin = ROOT.gStyle.GetPadBottomMargin()
    usable_width = 1. - l_margin - r_margin
    usable_height = 1. - t_margin - b_margin
    for i in xrange(len(split_points)+1):
        pad = ROOT.TPad('pad%i'%i, '', 0., 0., 1., 1.)
        pad_l_margin = l_margin + sum(split_points[0:i]) * usable_width + gaps_left[i-1]
        pad_r_margin = (1. - sum(split_points[0:i+1])) * usable_width + r_margin + gaps_right[i]
        if i > 0:
            pad.SetLeftMargin(pad_l_margin)
        if i < len(split_points):
            pad.SetRightMargin(pad_r_margin)
        print pad.GetLeftMargin(), pad.GetRightMargin()
        if ratio is not None:
            pad.SetBottomMargin(b_margin + usable_height * ratio)
        pad.SetFillStyle(4000)
        pad.Draw()
        pads.append(pad)
        if ratio is not None:
            padr = ROOT.TPad('r_pad%i'%i, '', 0., 0., 1., 1.)
            if i > 0:
                padr.SetLeftMargin(pad_l_margin)
            if i < len(split_points):
                padr.SetRightMargin(pad_r_margin)
            padr.SetTopMargin(1 - (b_margin + usable_height * ratio))
            padr.SetFillStyle(4000)
            padr.Draw()
            ratio_pads.append(padr)


    pads[0].cd()
    # pads.reverse()
    return pads, ratio_pads



years = args.years.split(',')
chg = args.charge

ref_hists = {}

# These factors undo the kNNLO that is applied for the nominal prediction
corr_factors = {
    '2016': 178.9 / 214.68,
    '2017': 191.4 / 229.92,
    '2018': 191.4 / 229.92
}
# output_2018_fid_region_fid_pt_binned.root

doScaleUncert = True

for yr in years:
    f_ref = ROOT.TFile('output_%s_%s_%s.root' % (yr, args.selection, args.scheme))
    h_ref_m = Node()
    TDirToNode(f_ref, 'm/XS/2D', h_ref_m)
    h_ref_e = Node()
    TDirToNode(f_ref, 'e/XS/2D', h_ref_e)

    if doScaleUncert:
        h_xs_m = DoScaleEnvelope(h_ref_m, 'XS_WG_%s_m_acc' % chg)
        h_xs_e = DoScaleEnvelope(h_ref_e, 'XS_WG_%s_e_acc' % chg)
    else:
        h_xs_m = h_ref_m['XS_WG_%s_m_acc' % chg]
        h_xs_e = h_ref_e['XS_WG_%s_e_acc' % chg]

    h_xs = h_xs_m.Clone()
    h_xs.Add(h_xs_e)

    # Scale from pb to fb, then a factor 3/2 since we only added e + mu, and we want all l
    h_xs.Scale(1. * (3. / 2.))
    # h_xs.Scale(1000. * (3. / 2.))
    h_xs.Scale(corr_factors[yr])
    h_xs.Print("range")
    h_xs.Scale(1, 'width')
    h_xs.Print("range")

    ref_hists[yr] = h_xs

n_bins_pt = ref_hists[years[0]].GetNbinsX()
n_bins_phi = ref_hists[years[0]].GetNbinsY()
width = 1. / n_bins_phi

canv = ROOT.TCanvas('diff_xsec', 'diff_xsec')

ratio = None
if args.ratio:
    ratio = 0.4
pads, ratio_pads = SetupPads([width] * (n_bins_phi - 1), [0, 0], [0, 0], ratio=ratio)

xsec_results = {}
with open('%s.json' % args.scheme) as jsonfile:
    xsec_results = json.load(jsonfile)['xsec2D']

ref_hists_1D = []
obs_graphs = []

for i in range(n_bins_phi):
    ref_hists_1D.append(ref_hists[yr].ProjectionX('proj_%i' % i, i + 1, i + 1, 'e'))
    h_1D = ref_hists_1D[-1]
    x_vals = []
    y_vals = []
    ex_hi = []
    ex_lo = []
    ey_hi = []
    ey_lo = []

    for j in range(n_bins_pt):
        x_vals.append(h_1D.GetXaxis().GetBinCenter(j + 1))
        # ex_lo.append(h_1D.GetXaxis().GetBinWidth(j + 1) / 2.)
        # ex_hi.append(h_1D.GetXaxis().GetBinWidth(j + 1) / 2.)
        ex_lo.append(0.)
        ex_hi.append(0.)
        poi_res = xsec_results['r_%s_%i' % (chg, j)]
        y_vals.append(h_1D.GetBinContent(j + 1) * poi_res['Val'])
        ey_hi.append(h_1D.GetBinContent(j + 1) * poi_res['ErrorHi'])
        ey_lo.append(h_1D.GetBinContent(j + 1) * poi_res['ErrorLo'] * -1.)
    obs_graphs.append(ROOT.TGraphAsymmErrors(len(x_vals), array('d', x_vals), array('d', y_vals), array('d', ex_lo), array('d', ex_hi), array('d', ey_lo), array('d', ey_hi)))

h_axes = [h.Clone() for h in ref_hists_1D]
r_h_axes = [h.Clone() for h in ref_hists_1D]

legend = ROOT.TLegend(0.74, 0.74, 0.94, 0.86, '', 'NBNDC')

h_obs = ROOT.TH1F('h_obs', '', 1, 0, 1)
plot.Set(h_obs, LineWidth=2)

latex = ROOT.TLatex()
latex.SetTextFont(62)
latex.SetTextSize(0.04)
latex.SetTextAlign(22)
# latex.SetTextColor(14)

for i, h in enumerate(h_axes):
    h.Reset()
    # h.SetMinimum(5E-5)
    # h.SetMaximum(8E-1)
    h.SetMinimum(1E-2)
    h.SetMaximum(5E-0)
    pads[i].cd()
    pads[i].SetLogy(True)
    h.Draw()
    h.GetXaxis().SetNdivisions(510)
    h.GetXaxis().ChangeLabel(-1, -1., -1., -1, -1, -1, ' ')
    h.GetYaxis().SetTickLength(h.GetYaxis().GetTickLength() * 0.5)
    if i > 0:
        h.GetYaxis().SetLabelSize(0)
    if i == n_bins_phi - 1:
        h.GetXaxis().SetTitle('Photon p_{T} (GeV)')
    if i == 0:
        h.GetYaxis().SetTitle('d^{2}#sigma/dp_{T}^{#gamma}d|#phi_{f}| (fb/GeV)')

    plot.Set(ref_hists_1D[i], LineColor=2, FillColor=plot.CreateTransparentColor(2, 0.2))
    plot.Set(obs_graphs[i], LineWidth=2)
    ref_hists_1D[i].Draw('E2SAME')
    ref_hists_1D[i].Draw('LSAME')
    obs_graphs[i].Draw('SAMEP')

    if i == 0:
        legend.AddEntry(h_obs, 'Observed', 'PE')
        legend.AddEntry(ref_hists_1D[i], 'W#gamma NLO (MG5_aMC@NLO)', 'LF')

    pad_width = 1. - pads[i].GetLeftMargin() - pads[i].GetRightMargin()
    latex.DrawLatexNDC(pads[i].GetLeftMargin() + pad_width / 2., 0.9, '%.2f #leq |#phi_{f}| < %.2f' % (ref_hists[yr].GetYaxis().GetBinLowEdge(i + 1), ref_hists[yr].GetYaxis().GetBinUpEdge(i + 1)))

r_ref_hists_1D = []
r_obs_graphs = []

if args.ratio:
    for i, h in enumerate(r_h_axes):
        h.Reset()
        # h.SetMinimum(5E-5)
        # h.SetMaximum(8E-1)
        h.SetMinimum(0.8)
        h.SetMaximum(1.2)
        ratio_pads[i].cd()
        h.GetXaxis().SetNdivisions(510)
        h.GetXaxis().ChangeLabel(-1, -1., -1., -1, -1, -1, ' ')
        h.Draw()
        r_ref_hists_1D.append(plot.MakeRatioHist(ref_hists_1D[i], ref_hists_1D[i], True, False))
        r_obs_graphs.append(DivideGraphByHist(obs_graphs[i], ref_hists_1D[i]))
        r_ref_hists_1D[i].Draw('E2SAME')
        r_ref_hists_1D[i].SetMarkerSize(0)
        r_obs_graphs[i].Draw('SAMEP')
        # r_ref_hists_1D[i].Draw('LSAME')

latex.SetTextSize(0.05)
chg_labels = {
    'p': '+',
    'n': '-',
    'x': '#pm'
}
latex.DrawLatexNDC(0.17, 0.17, 'W^{%s}(#rightarrowl^{%s}#nu)#gamma' % (chg_labels[args.charge], chg_labels[args.charge]))

legend.Draw()

# pads[0].cd()
# pads[1].cd()
# ref_hists[yr].ProjectionX('proj_1', 2, 2, 'e').Draw()
# pads[2].cd()
# ref_hists[yr].ProjectionX('proj_2', 3, 3, 'e').Draw()

# pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

# h_xs = ref_hists[years[0]]

# h_axes = [h_xs.Clone() for x in pads]
# for h in h_axes:
#     h.Reset()
# h_axes[1].GetXaxis().SetTitle('p_{T}^{#gamma} (GeV)')
# h_axes[0].GetYaxis().SetTitle('#Delta#sigma/#Deltap_{T}^{#gamma} (fb/GeV)')
# pads[0].SetLogy()
# h_axes[0].SetMinimum(1E-5)
# h_axes[0].Draw()

# plot.Set(ref_hists['2016'], LineWidth=1, LineColor=8, MarkerColor=8)
# plot.Set(ref_hists['2017'], LineWidth=1, LineColor=2, MarkerColor=2)
# plot.Set(ref_hists['2018'], LineWidth=1, LineColor=4, MarkerColor=4)

# ref_hists['2016'].Draw('HISTSAMEE')
# ref_hists['2017'].Draw('HISTSAMEE')
# ref_hists['2018'].Draw('HISTSAMEE')

# f_fit = ROOT.TFile('multidimfit.root')
# fitres = f_fit.Get('fit_mdf')

# h_res = ref_hists['2018'].Clone()
# for ib in xrange(1, h_res.GetNbinsX() + 1):
#     bin_xs = h_res.GetBinContent(ib)
#     var = fitres.randomizePars().find('r_%s_%i' % (chg, ib - 1))
#     h_res.SetBinContent(ib, var.getVal() * bin_xs)
#     h_res.SetBinError(ib, var.getError() * bin_xs)

# plot.Set(h_res, LineWidth=2, LineColor=1, MarkerColor=1)

# h_res.Draw('HISTSAMEE')
# h_axes[0].SetMinimum(1E-7)

# plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.43)

# pads[1].cd()
# pads[1].SetGrid(0, 1)
# h_axes[1].Draw()

# r_xs_2016 = plot.MakeRatioHist(ref_hists['2016'], ref_hists['2018'], True, False)
# r_xs_2017 = plot.MakeRatioHist(ref_hists['2017'], ref_hists['2018'], True, False)
# r_xs_2018 = plot.MakeRatioHist(ref_hists['2018'], ref_hists['2018'], True, False)
# r_res = plot.MakeRatioHist(h_res, ref_hists['2018'], True, False)

# r_xs_2018.SetFillColor(plot.CreateTransparentColor(12, 0.3))
# r_xs_2018.SetMarkerSize(0)
# r_xs_2018.Draw('E2SAME')
# r_xs_2016.Draw('SAMEE')
# r_xs_2017.Draw('SAMEE')
# r_res.Draw('SAMEE')

# plot.SetupTwoPadSplitAsRatio(
#     pads, plot.GetAxisHist(
#         pads[0]), plot.GetAxisHist(pads[1]), 'Ratio', True, 0.61, 1.39)

# pads[0].cd()
# pads[0].GetFrame().Draw()
# pads[0].RedrawAxis()

pads[0].cd()
plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.24, 0.035, 1.2, cmsTextSize=0.9)
plot.DrawTitle(pads[0], '136.9 fb^{-1} (13 TeV)', 3)


canv.Print('.png')
canv.Print('.pdf')
