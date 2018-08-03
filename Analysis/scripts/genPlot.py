import ROOT
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

tname = 'WGAnalysis'

samples = {
    'WG': 'output/100718/wgamma_2016_v1/wg_gen/WGToLNuG-EFT_pTA_300_inf-madgraphMLM_0.root'
}

remap = {
    'WG': 'WGToLNuG-EFT_pTA_300_inf-madgraphMLM'
}


hists = Node()

for name, wt in [
        ('nominal', '1.0'),
        ('C3w_0p1', 'weight_C3w_0p1'),
        ('C3w_0p2', 'weight_C3w_0p2'),
        ('C3w_0p4', 'weight_C3w_0p4'),
        ('C3w_0p8', 'weight_C3w_0p8'),
        ('C3w_1p6', 'weight_C3w_1p6')]:
    hists[name] = Hist('TH1D', (10, -3.15, 3.15), 'WG', ['phi1'],
        sel='g_pt>300. && l_pt>80. && n_pt>80. && fabs(l_eta) < 2.4 && l_g_dr > 3.0 && fabs(g_eta) < 3. && weight_C3w_1p6 < 10.', wt=wt)

MultiDraw(hists, samples, tname)

hists.ForEach(lambda x: NormaliseTo(x, 1.0))

canv = ROOT.TCanvas('gen_plot', 'gen_plot')
pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

# Get the data and create axis hist
h_nominal = hists['nominal']
h_0p2 = hists['C3w_0p2']
h_0p4 = hists['C3w_0p4']
h_0p8 = hists['C3w_0p8']
h_1p6 = hists['C3w_1p6']

h_axes = [h_nominal.Clone() for x in pads]
for h in h_axes:
    # h.GetXaxis().SetLimits(2.1,200)
    h.Reset()

h_axes[1].GetXaxis().SetTitle('#varphi')
h_axes[0].GetYaxis().SetTitle('a.u.')
h_axes[0].Draw()

# A dict to keep track of the hists
legend = ROOT.TLegend(0.67, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

legend.AddEntry(h_nominal, 'Nominal', 'L')
legend.AddEntry(h_0p2, 'C_{3W} = 0.2', 'L')
# legend.AddEntry(h_0p4, 'C_{3W} = 0.4', 'L')
legend.AddEntry(h_0p8, 'C_{3W} = 0.8', 'L')
legend.AddEntry(h_1p6, 'C_{3W} = 1.6', 'L')

plot.Set(h_0p2, LineColor=2, LineWidth=2)
plot.Set(h_0p4, LineColor=4, LineWidth=2)
plot.Set(h_0p8, LineColor=8, LineWidth=2)
plot.Set(h_1p6, LineColor=9, LineWidth=2)

h_nominal.Draw('HISTSAME')
h_0p2.Draw('HISTSAME')
# h_0p4.Draw('HISTSAME')
h_0p8.Draw('HISTSAME')
h_1p6.Draw('HISTSAME')

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.35)
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
r_0p2.Draw('HISTSAME')
# r_0p2.Draw('HISTSAME')
r_0p8.Draw('HISTSAME')
r_1p6.Draw('HISTSAME')
# r_data.Draw('SAME')
plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'C_{3W} = X / Nominal', True, 0.81, 1.19)


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
# plot.DrawTitle(pads[0], args.title, 1)

# ... and we're done
canv.Print('.png')
canv.Print('.pdf')

fout = ROOT.TFile('output_gen.root', 'RECREATE')

for path, name, obj in hists.ListObjects():
    print (path, name, obj)
    WriteToTFile(obj, fout, path, name)

fout.Close()
