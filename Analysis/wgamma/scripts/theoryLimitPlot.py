import ROOT
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

# parser = argparse.ArgumentParser()
# parser.add_argument('--draw', default='phi1')
# parser.add_argument('--abs', action='store_true')
# parser.add_argument('--unit-norm', action='store_true')
# parser.add_argument('--charge', default='+1')
# parser.add_argument('--g_pt', default='300.')
# parser.add_argument('--g_eta', default='3.')
# parser.add_argument('--l_pt', default='80.')
# parser.add_argument('--l_eta', default='2.4')
# parser.add_argument('--n_pt', default='80.')
# parser.add_argument('--n_eta', default='9999.')
# parser.add_argument('--dr', default='3.0')
# parser.add_argument('--output', '-o', default='gen_plot')
# # parser.add_argument('--ratio', '-o', default='gen_plot')

# args = parser.parse_args()

def ConvertToHist(gr, bins):
    N = gr.GetN()
    hist = ROOT.TH1F('', '', N, array('d', pt_bins[0:N + 1]))
    for i in xrange(N):
        hist.SetBinContent(i + 1, gr.GetY()[i])
    gr.Print()
    hist.Print("range")
    return hist

canv = ROOT.TCanvas('theory_limits', 'theory_limits')
pads = plot.OnePad()

pt_bins = [150, 210, 300, 420, 600, 850, 1200]

gr_0jet = plot.LimitTGraphFromJSONFile('limit_inc_bsm.json', 'obs')
h_0jet = ConvertToHist(gr_0jet, pt_bins)

gr_0jet_1bin = plot.LimitTGraphFromJSONFile('limit_inc_1bin_bsm.json', 'obs')
h_0jet_1bin = ConvertToHist(gr_0jet_1bin, pt_bins)

gr_inc = plot.LimitTGraphFromJSONFile('limit_inc.json', 'obs')
h_inc = ConvertToHist(gr_inc, pt_bins)

gr_inc_1bin = plot.LimitTGraphFromJSONFile('limit_inc_1bin.json', 'obs')
h_inc_1bin = ConvertToHist(gr_inc_1bin, pt_bins)

pads[0].SetLogy(1)
h_axes = [h_0jet.Clone() for x in pads]
for h in h_axes:
    h.Reset()

plot.Set(h_axes[0], Minimum=0.01, Maximum=50)
plot.Set(h_axes[0].GetXaxis(), Title='Maximum p_{T}^{#gamma} (GeV)')
plot.Set(h_axes[0].GetYaxis(), Title='95% CL limit on C_{3W} (TeV^{-2})')
h_axes[0].Draw()

plot.Set(h_0jet, LineColor=1, LineWidth=2)
plot.Set(h_0jet_1bin, LineColor=1, LineWidth=2, LineStyle=7)
plot.Set(h_inc, LineColor=2, LineWidth=2)
plot.Set(h_inc_1bin, LineColor=2, LineWidth=2, LineStyle=7)

h_0jet.Draw('SAME')
h_0jet_1bin.Draw('SAME')
h_inc.Draw('SAME')
h_inc_1bin.Draw('SAME')

legend = ROOT.TLegend(0.45, 0.7, 0.90, 0.91, '', 'NBNDC')
legend.AddEntry(h_0jet, 'LO, #leq3 jets', 'L')
legend.AddEntry(h_0jet_1bin, 'LO, #leq3 jets, no #varphi binning', 'L')
legend.AddEntry(h_inc, 'LO, #leq3 jets, No BSM', 'L')
legend.AddEntry(h_inc_1bin, 'LO, #leq3 jets, No BSM, no #varphi binning', 'L')
legend.Draw()

plot.DrawTitle(pads[0], '13 TeV', 3)
plot.DrawTitle(pads[0], '136.7 fb^{-1}: W^{+}(e^{+}#nu/#mu^{+}#nu)#gamma', 1)

# plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)

# Get the data and create axis hist
# h_nominal = hists['nominal']
# h_sm = hists['SM']
# h_th = hists['TH']
# h_0p1 = hists['C3w_0p1']
# h_0p2 = hists['C3w_0p2']
# h_0p4 = hists['C3w_0p4']
# h_1p0 = hists['C3w_1p0']




# if 'lhe_phi1' in drawvar:
#     labelvar = '#varphi_{true}'
# else:
#     labelvar = '#varphi_{reco}'
# if args.abs:
#     labelvar = '|%s|' % labelvar
# h_axes[1].GetXaxis().SetTitle(labelvar)

# h_axes[0].GetYaxis().SetTitle('a.u.')
# h_axes[0].Draw()

# # A dict to keep track of the hists
# legend = ROOT.TLegend(0.67, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

# legend.AddEntry(h_nominal, 'EWDim6', 'L')
# legend.AddEntry(h_sm, 'SM', 'L')
# legend.AddEntry(h_th, 'Reference', 'L')
# legend.AddEntry(h_0p1, 'C_{3W} = 0.1', 'L')
# legend.AddEntry(h_0p2, 'C_{3W} = 0.2', 'L')
# legend.AddEntry(h_0p4, 'C_{3W} = 0.4', 'L')
# legend.AddEntry(h_1p0, 'C_{3W} = 1.0', 'L')

# plot.Set(h_sm, LineColor=2, LineWidth=2, MarkerColor=2)
# plot.Set(h_th, LineColor=4, LineWidth=2, MarkerColor=4)
# plot.Set(h_0p1, LineColor=6, LineWidth=1)
# plot.Set(h_0p2, LineColor=7, LineWidth=1)
# plot.Set(h_0p4, LineColor=9, LineWidth=1)
# plot.Set(h_1p0, LineColor=28, LineWidth=1)

# h_nominal.Draw('HISTSAMEE')
# # h_sm.Draw('HISTSAMEE')
# # h_th.Draw('HISTSAMEE')
# h_0p1.Draw('HISTSAME')
# h_0p2.Draw('HISTSAME')
# h_0p4.Draw('HISTSAME')
# h_1p0.Draw('HISTSAME')

# plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.43)
# legend.Draw()
# # plot.FixBoxPadding(pads[0], legend, 0.05)

# # # Do the ratio plot
# pads[1].cd()
# pads[1].SetGrid(0, 1)
# h_axes[1].Draw()

# # r_data = plot.MakeRatioHist(h_data, h_tot, True, False)
# # r_nominal = plot.MakeRatioHist(h_nominal, h_nominal, True, False)
# r_0p1 = plot.MakeRatioHist(h_0p1, h_nominal, True, False)
# r_0p2 = plot.MakeRatioHist(h_0p2, h_nominal, True, False)
# r_0p4 = plot.MakeRatioHist(h_0p4, h_nominal, True, False)
# r_1p0 = plot.MakeRatioHist(h_1p0, h_nominal, True, False)
# r_0p1.Draw('HISTSAME')
# r_0p2.Draw('HISTSAME')
# r_0p4.Draw('HISTSAME')
# r_1p0.Draw('HISTSAME')
# # r_data.Draw('SAME')
# plot.SetupTwoPadSplitAsRatio(
#     pads, plot.GetAxisHist(
#         pads[0]), plot.GetAxisHist(pads[1]), 'C_{3W} = X / Nominal', True, 0.64, 1.36)


# # Go back and tidy up the axes and frame
# pads[0].cd()
# pads[0].GetFrame().Draw()
# pads[0].RedrawAxis()

# # CMS logo
# # plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)
# # plot.DrawTitle(pads[0], '0.1 fb^{-1} (13 TeV)', 3)

# latex = ROOT.TLatex()
# plot.Set(latex, NDC=None, TextFont=42, TextSize=0.03)
# # latex.DrawLatex(0.20, 0.75, args.title)
# plot.DrawTitle(pads[0], '#sqrt{s} = 13 TeV', 3)
# plot.DrawTitle(pads[0], 'W^{%s}#gamma^{} LO - MG5_aMC@NLO 2.4.2' %
#                ('+' if args.charge == '+1' else '-' if args.charge == '-1' else '#pm'), 1)

# pt_l = ROOT.TPaveText(0.23, 0.65, 0.55, 0.9, 'NDCNB')
# pt_l.AddText('Selection: before showering (LHE)')
# pt_l.AddText('p_{T}^{#gamma} > %s GeV, |#eta^{#gamma}| < %s' %
#              (args.g_pt, args.g_eta))
# pt_l.AddText('p_{T}^{l} > %s GeV, |#eta^{l}| < %s' % (args.l_pt, args.l_eta))
# pt_l.AddText('p_{T}^{miss} > %s GeV' % args.n_pt)
# pt_l.AddText('#DeltaR(l, #gamma) > %s' % args.dr)
# pt_l.AddText('nJets (ME) == 0')
# plot.Set(pt_l, TextAlign=11, TextFont=42, BorderSize=0)

# pt_l.Draw()
# ... and we're done
canv.Print('.png')
canv.Print('.pdf')

# fout = ROOT.TFile('output_gen.root', 'RECREATE')

# fout.Close()
