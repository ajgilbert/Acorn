import ROOT
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

tname = 'WGAnalysis'

samples = {
    'WGSM': 'output/100718/wgamma_2016_v2/wg_gen/WpAToLNuA0j_5f_LO_MLM_pTA_300_0.root',
    'WGTH': 'output/100718/wgamma_2016_v2/wg_gen/Theorists_0.root',
    'WG': 'output/100718/wgamma_2016_v2/wg_gen/WpAToLNuA0j_5f_LO_MLM_EFT_pTA_300_0.root'
}

remap = {
    'WG': 'WGToLNuG-EFT_pTA_300_inf-madgraphMLM'
}


parser = argparse.ArgumentParser()
parser.add_argument('--draw', default='phi1')
parser.add_argument('--abs', action='store_true')
parser.add_argument('--charge', default='+1')
parser.add_argument('--g_pt', default='300.')
parser.add_argument('--g_eta', default='3.')
parser.add_argument('--l_pt', default='80.')
parser.add_argument('--l_eta', default='2.4')
parser.add_argument('--n_pt', default='80.')
parser.add_argument('--n_eta', default='9999.')
parser.add_argument('--dr', default='3.0')
parser.add_argument('--output', '-o', default='gen_plot')
# parser.add_argument('--ratio', '-o', default='gen_plot')

args = parser.parse_args()


hists = Node()

# drawabs = True
# drawvar = 'lhe_phi1'
if args.abs:
    drawvar = 'fabs(%s)' % args.draw
    binning = (5, 0, 3.15)
else:
    drawvar = '%s' % args.draw
    binning = (10, -3.15, 3.15)


for name, sa, wt in [
        ('nominal', 'WG', '1.0'),
        ('C3w_0p1', 'WG', 'weight_C3w_0p1'),
        ('C3w_0p2', 'WG', 'weight_C3w_0p2'),
        ('C3w_0p4', 'WG', 'weight_C3w_0p4'),
        ('C3w_0p8', 'WG', 'weight_C3w_0p8'),
        ('C3w_1p6', 'WG', 'weight_C3w_1p6'),
        ('SM',      'WGSM', '1.0'),
        ('TH',      'WGTH', '1.0')
        ]:
    # hists[name] = Hist('TH1D', (10, -3.15, 3.15), 'WG', [drawvar],
    sel = 'l_charge==%s && nparts>=3 && nparts <=3 && g_pt>%s && l_pt>%s && n_pt>%s && fabs(l_eta) < %s && fabs(n_eta) < %s && fabs(g_eta) < %s && l_g_dr > %s && weight_C3w_0p2 < 10.' % (args.charge, args.g_pt, args.l_pt, args.n_pt, args.l_eta, args.n_eta, args.g_eta, args.dr)
    print sel
    hists[name] = Hist('TH1D', binning, sa, [drawvar],
        # sel='l_charge==+1 && nparts>=3 && nparts <=3 && g_pt>300. && l_pt>80. && n_pt>80. && fabs(l_eta) < 2.4 && l_g_dr > 3.0 && fabs(g_eta) < 3.', wt=wt)
        sel=sel, wt=wt)


MultiDraw(hists, samples, tname)

hists.ForEach(lambda x: NormaliseTo(x, 1.0))
# hists.ForEach(lambda x: WidthDivide(x))

canv = ROOT.TCanvas(args.output, args.output)
pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

# Get the data and create axis hist
h_nominal = hists['nominal']
h_sm = hists['SM']
h_th = hists['TH']
h_0p2 = hists['C3w_0p2']
h_0p4 = hists['C3w_0p4']
h_0p8 = hists['C3w_0p8']
h_1p6 = hists['C3w_1p6']

h_axes = [h_nominal.Clone() for x in pads]
for h in h_axes:
    # h.GetXaxis().SetLimits(2.1,200)
    h.Reset()


if 'lhe_phi1' in drawvar:
    labelvar = '#varphi_{true}'
else:
    labelvar = '#varphi_{reco}'
if args.abs:
    labelvar = '|%s|' % labelvar
h_axes[1].GetXaxis().SetTitle(labelvar)

h_axes[0].GetYaxis().SetTitle('a.u.')
h_axes[0].Draw()

# A dict to keep track of the hists
legend = ROOT.TLegend(0.67, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

legend.AddEntry(h_nominal, 'EWDim6', 'L')
legend.AddEntry(h_sm, 'SM', 'L')
legend.AddEntry(h_th, 'Reference', 'L')
# legend.AddEntry(h_0p2, 'C_{3W} = 0.2', 'L')
# legend.AddEntry(h_0p4, 'C_{3W} = 0.4', 'L')
# legend.AddEntry(h_0p8, 'C_{3W} = 0.8', 'L')
# legend.AddEntry(h_1p6, 'C_{3W} = 1.6', 'L')

plot.Set(h_sm, LineColor=2, LineWidth=2, MarkerColor=2)
plot.Set(h_th, LineColor=4, LineWidth=2, MarkerColor=4)
plot.Set(h_0p2, LineColor=8, LineWidth=1)
plot.Set(h_0p4, LineColor=4, LineWidth=2)
plot.Set(h_0p8, LineColor=8, LineWidth=2)
plot.Set(h_1p6, LineColor=9, LineWidth=2)

h_nominal.Draw('HISTSAMEE')
h_sm.Draw('HISTSAMEE')
h_th.Draw('HISTSAMEE')
# h_0p2.Draw('HISTSAME')
# h_0p4.Draw('HISTSAME')
# h_0p8.Draw('HISTSAME')
# h_1p6.Draw('HISTSAME')

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.43)
legend.Draw()
# plot.FixBoxPadding(pads[0], legend, 0.05)

# # Do the ratio plot
pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].Draw()

# r_data = plot.MakeRatioHist(h_data, h_tot, True, False)
# r_nominal = plot.MakeRatioHist(h_nominal, h_nominal, True, False)
r_0p2 = plot.MakeRatioHist(h_0p2, h_nominal, True, False)
r_0p4 = plot.MakeRatioHist(h_0p4, h_nominal, True, False)
r_0p8 = plot.MakeRatioHist(h_0p8, h_nominal, True, False)
r_1p6 = plot.MakeRatioHist(h_1p6, h_nominal, True, False)
# r_0p2.Draw('HISTSAME')
# r_0p2.Draw('HISTSAME')
# r_0p8.Draw('HISTSAME')
# r_1p6.Draw('HISTSAME')
# r_data.Draw('SAME')
plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'C_{3W} = X / Nominal', True, 0.74, 1.26)


# Go back and tidy up the axes and frame
pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

# CMS logo
# plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)
# plot.DrawTitle(pads[0], '0.1 fb^{-1} (13 TeV)', 3)

latex = ROOT.TLatex()
plot.Set(latex, NDC=None, TextFont=42, TextSize=0.03)
# latex.DrawLatex(0.20, 0.75, args.title)
plot.DrawTitle(pads[0], '#sqrt{s} = 13 TeV', 3)
plot.DrawTitle(pads[0], 'W^{%s}#gamma^{} LO - MG5_aMC@NLO 2.4.2' % ('+' if args.charge == '+1' else '-'), 1)

pt_l = ROOT.TPaveText(0.23, 0.65, 0.55, 0.9, 'NDCNB')
pt_l.AddText('Selection: before showering (LHE)')
pt_l.AddText('p_{T}^{#gamma} > %s GeV, |#eta^{#gamma}| < %s' % (args.g_pt, args.g_eta))
pt_l.AddText('p_{T}^{l} > %s GeV, |#eta^{l}| < %s' % (args.l_pt, args.l_eta))
pt_l.AddText('p_{T}^{miss} > %s GeV' % args.n_pt)
pt_l.AddText('#DeltaR(l, #gamma) > %s' % args.dr)
pt_l.AddText('nJets (ME) == 0')
plot.Set(pt_l, TextAlign=11, TextFont=42, BorderSize=0)

pt_l.Draw()
# ... and we're done
canv.Print('.png')
canv.Print('.pdf')

fout = ROOT.TFile('output_gen.root', 'RECREATE')

for path, name, obj in hists.ListObjects():
    print (path, name, obj)
    WriteToTFile(obj, fout, path, name)

fout.Close()
