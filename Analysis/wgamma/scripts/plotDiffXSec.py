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
plot.ModTDRStyle()

years = ['2016', '2017', '2018']
chg = 'p'

ref_hists = {}

# These factors undo the kNNLO that is applied for the nominal prediction
corr_factors = {
    '2016': 178.9 / 214.68,
    '2017': 191.4 / 229.92,
    '2018': 191.4 / 229.92
}

for yr in years:
    f_ref = ROOT.TFile('output_%s_eft_region_default.root' % yr)
    h_ref_m = Node()
    TDirToNode(f_ref, 'm/XS/gen_p0_pt', h_ref_m)
    h_ref_e = Node()
    TDirToNode(f_ref, 'e/XS/gen_p0_pt', h_ref_e)

    h_xs_m = h_ref_m['XS_WG_%s_m_acc' % chg]
    h_xs_e = h_ref_e['XS_WG_%s_e_acc' % chg]

    h_xs = h_xs_m.Clone()
    h_xs.Add(h_xs_e)

    # Scale from pb to fb, then a factor 3/2 since we only added e + mu, and we want all l
    h_xs.Scale(1000. * (3. / 2.))
    h_xs.Scale(corr_factors[yr])
    h_xs.Scale(1, 'width')

    ref_hists[yr] = h_xs

canv = ROOT.TCanvas('diff_xsec', 'diff_xsec')
pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

h_xs = ref_hists[years[0]]

h_axes = [h_xs.Clone() for x in pads]
for h in h_axes:
    h.Reset()
h_axes[1].GetXaxis().SetTitle('p_{T}^{#gamma} (GeV)')
h_axes[0].GetYaxis().SetTitle('#Delta#sigma/#Deltap_{T}^{#gamma} (fb/GeV)')
pads[0].SetLogy()
h_axes[0].SetMinimum(1E-5)
h_axes[0].Draw()

plot.Set(ref_hists['2016'], LineWidth=1, LineColor=8, MarkerColor=8)
plot.Set(ref_hists['2017'], LineWidth=1, LineColor=2, MarkerColor=2)
plot.Set(ref_hists['2018'], LineWidth=1, LineColor=4, MarkerColor=4)

ref_hists['2016'].Draw('HISTSAMEE')
ref_hists['2017'].Draw('HISTSAMEE')
ref_hists['2018'].Draw('HISTSAMEE')

f_fit = ROOT.TFile('multidimfit.root')
fitres = f_fit.Get('fit_mdf')

h_res = ref_hists['2018'].Clone()
for ib in xrange(1, h_res.GetNbinsX() + 1):
    bin_xs = h_res.GetBinContent(ib)
    var = fitres.randomizePars().find('r_%s_%i' % (chg, ib - 1))
    h_res.SetBinContent(ib, var.getVal() * bin_xs)
    h_res.SetBinError(ib, var.getError() * bin_xs)

plot.Set(h_res, LineWidth=2, LineColor=1, MarkerColor=1)

h_res.Draw('HISTSAMEE')
h_axes[0].SetMinimum(1E-7)

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.43)

pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].Draw()

r_xs_2016 = plot.MakeRatioHist(ref_hists['2016'], ref_hists['2018'], True, False)
r_xs_2017 = plot.MakeRatioHist(ref_hists['2017'], ref_hists['2018'], True, False)
r_xs_2018 = plot.MakeRatioHist(ref_hists['2018'], ref_hists['2018'], True, False)
r_res = plot.MakeRatioHist(h_res, ref_hists['2018'], True, False)

r_xs_2018.SetFillColor(plot.CreateTransparentColor(12, 0.3))
r_xs_2018.SetMarkerSize(0)
r_xs_2018.Draw('E2SAME')
r_xs_2016.Draw('SAMEE')
r_xs_2017.Draw('SAMEE')
r_res.Draw('SAMEE')

plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'Ratio', True, 0.61, 1.39)

pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

canv.Print('.png')
canv.Print('.pdf')
