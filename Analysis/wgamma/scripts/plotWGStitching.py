import ROOT
# import json
# import math
from pprint import pprint
from collections import defaultdict
from Acorn.Analysis.analysis import *
import CombineHarvester.CombineTools.plotting as plot

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

f_inc = ROOT.TFile('dist_lhe_all_nofilter_total.root')
h_inc = f_inc.Get('h')
n_pos_inc = 296482
n_neg_inc = 93518
xsec_inc = 3.080e+02

h_inc.Scale(1. / (n_pos_inc + n_neg_inc))

f_130 = ROOT.TFile('dist_lhe_pt130_nofilter.root')
h_130 = f_130.Get('h')
n_pos_130 = 15230
n_neg_130 = 4770
xsec_130 = 2.1392269
h_130.Scale(1. / (n_pos_130 + n_neg_130))

canv = ROOT.TCanvas('wg_stitching', 'wg_stitching')
pads = plot.TwoPadSplit(0.35, 0.01, 0.01)

h_axes = [h_130.Clone() for x in pads]
for h in h_axes:
    h.Reset()

h_axes[0].GetYaxis().SetTitle('a.u.')
h_axes[0].Draw()
# pads[0].SetLogy()
# # pads[0].SetLogx()
# # pads[1].SetLogx()
h_axes[0].SetMinimum(-0.2)
h_axes[0].SetMaximum(2.0)
h_axes[0].GetXaxis().SetRangeUser(50, 200)
h_axes[0].GetXaxis().SetTitle('LHE Photon p_{T} [GeV]')

# # # # A dict to keep track of the hists
legend = ROOT.TLegend(0.6, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

# h_pdf_env_fill = h_pdf_env.Clone()
plot.Set(h_130, LineColor=2, LineWidth=2, MarkerColor=2)
plot.Set(h_inc, LineColor=1, LineWidth=1, MarkerColor=1)
# plot.Set(h_pdf_env_fill, LineColor=2, LineWidth=2, MarkerColor=2, MarkerSize=0, FillColorAlpha=(2, 0.3))

legend.AddEntry(h_130, 'p_{T}^{#gamma} > 120 GeV Gridpack', 'LF')
legend.AddEntry(h_inc, 'Inclusive Gridpack', 'L')

h_130.Draw('SAME')
h_130.Print('range')
h_inc.Draw('SAME')
h_inc.Print('range')
# # h_matrix_fill = h_matrix.Clone()
# # plot.Set(h_matrix, LineColor=4, LineWidth=2, MarkerColor=4)
# # plot.Set(h_matrix_fill, LineColor=4, LineWidth=2, MarkerColor=4, MarkerSize=0, FillColorAlpha=(4, 0.3))

# h_pdf_env_fill.Draw('E2SAME')
# # h_matrix_fill.Draw('E2SAME')
# h_pdf_env.Draw('HISTSAME')
# h_pdf_nlo.Draw('HISTSAME')
# # h_matrix.Draw('HISTSAME')
# h_axes[0].SetMinimum(1E-6)
# h_axes[0].SetMaximum(1E0)

# # plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)
legend.Draw()
# # # # plot.FixBoxPadding(pads[0], legend, 0.05)

# # # # # Do the ratio plot
pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].GetXaxis().SetRangeUser(50, 200)
h_axes[1].Draw()

r_MC = plot.MakeRatioHist(h_130, h_inc, True, True)
r_MC.Draw('SAME')
# r_pdf_env_fill = plot.MakeRatioHist(h_pdf_env_fill, h_pdf_env_fill, True, False)
# r_pdf_nlo = plot.MakeRatioHist(h_pdf_nlo, h_pdf_env_fill, True, False)
# r_pdf_env_fill.Draw('E2SAME')
# # r_matrix_fill.Draw('E2SAME')
# r_pdf_nlo.Draw('HISTSAME')
plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'Ratio', True, 0.5, 1.5)


# # # # Go back and tidy up the axes and frame
pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()


canv.Print('.png')
canv.Print('.pdf')
# fout = ROOT.TFile('output_pdf_uncertainty.root', 'RECREATE')

# NodeToTDir(fout, hists)

# fout.Close()
