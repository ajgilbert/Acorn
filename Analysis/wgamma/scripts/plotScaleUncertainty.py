import ROOT
# import json
import math
from pprint import pprint
from collections import defaultdict
import argparse
from Acorn.Analysis.analysis import *
import CombineHarvester.CombineTools.plotting as plot

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

def DoScaleEnvelope(node, nominal):
    h = node[nominal].Clone()
    h_alts = []
    for i in range(6):
        h_alts.append(node[nominal.replace('_nominal', '_scale_%i' % i)])
    for ix in xrange(1, h.GetNbinsX() + 1):
        for iy in xrange(1, h.GetNbinsY() + 1):
            max_dev = max([abs(x.GetBinContent(ix, iy) - h.GetBinContent(ix, iy)) for x in h_alts])
            h.SetBinError(ix, iy, max_dev)
    return h

parser = argparse.ArgumentParser()

parser.add_argument('--task', default='baseline', choices=['eft_region', 'baseline', 'photon_fakes'])
parser.add_argument('--indir', default='root://eoscms.cern.ch//store/cmst3/user/agilbert/191003-full/wgamma_2018_v4/WGamma_')

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

remap = {
    # 'WG-LO': 'WGToLNuG-madgraphMLM-stitched',
    # 'WG-LO-130': 'WGToLNuG-madgraphMLM-PtG-130',
    # 'WG-LO-500': 'WGToLNuG-madgraphMLM-PtG-500',
    'WG': 'WGToLNuG-amcatnloFXFX-stitched'
}

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')


hists = Node()
do_cats = []

# X = SelectionManager()

pt_binning = [30,40,50,60,80,100,150,200,300,500,800,1200]

drawvars = [
    ('gen_p0_pt', pt_binning),
]

npdf = 100

baseline_sel = 'gen_l0_pt>30 && abs(gen_l0_eta)<2.5 && gen_met>0 && gen_p0_pt>30 && abs(gen_p0_eta)<2.5 && gen_l0p0_dr>0.7 && lhe_frixione && gen_met>0'
for var, binning in drawvars:
    for sample in samples:
        hists['inclusive'][var][sample + '_nominal'] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel='%s && (gen_pdgid == 11 || gen_pdgid == 13)' % baseline_sel, wt='wt_def')
        for i in xrange(6):
            hists['inclusive'][var][sample + '_scale_%i' % i] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel='%s && (gen_pdgid == 11 || gen_pdgid == 13)' % baseline_sel, wt='wt_def * wt_sc_%i' % i)

MultiDraw(hists, samples, tname, mt_cores=4)

for var, binning in drawvars:
    for sample in samples:
        hists['inclusive'][var][sample+'_scale_env'] = DoScaleEnvelope(hists['inclusive'][var], sample + '_nominal')
        hists['inclusive'][var][sample+'_scale_env'].sample = sample
xsecs = {
    'WG': 178.9
}

nevents = {}

for sample in samples:
    f = ROOT.TFile.Open(samples[sample])
    nevents[sample] = f.Get('counters').GetBinContent(2)
    f.Close()


for path, hname, obj in hists.ListObjects():
    obj.Scale(1. * xsecs[obj.sample] / nevents[obj.sample])
    # We picked just one lepton flavour but will compare to all three
    # obj.Scale(3.)
    obj.Scale(1., 'width')


h_nominal = hists['inclusive']['gen_p0_pt']['WG_nominal'].Clone()
h_env = hists['inclusive']['gen_p0_pt']['WG_scale_env'].Clone()
h_scale_0 = hists['inclusive']['gen_p0_pt']['WG_scale_0'].Clone()
h_scale_1 = hists['inclusive']['gen_p0_pt']['WG_scale_1'].Clone()
h_scale_2 = hists['inclusive']['gen_p0_pt']['WG_scale_2'].Clone()
h_scale_3 = hists['inclusive']['gen_p0_pt']['WG_scale_3'].Clone()
h_scale_4 = hists['inclusive']['gen_p0_pt']['WG_scale_4'].Clone()
h_scale_5 = hists['inclusive']['gen_p0_pt']['WG_scale_5'].Clone()


canv = ROOT.TCanvas('scale_uncertainty', 'scale_uncertainty')
pads = plot.TwoPadSplit(0.35, 0.01, 0.01)

h_axes = [h_nominal.Clone() for x in pads]
for h in h_axes:
    h.Reset()

h_axes[0].GetYaxis().SetTitle('a.u.')
h_axes[0].Draw()
pads[0].SetLogy()
# pads[0].SetLogx()
# pads[1].SetLogx()
h_axes[0].SetMinimum(1E-9)
h_axes[0].GetXaxis().SetTitle('Photon p_{T} [GeV]')

# # # A dict to keep track of the hists
legend = ROOT.TLegend(0.67, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

plot.Set(h_nominal, LineColor=2, LineWidth=2, MarkerColor=2)
plot.Set(h_env, LineColor=2, LineWidth=2, MarkerColor=2, MarkerSize=0, FillColorAlpha=(2, 0.2))
plot.Set(h_scale_0, LineColor=8, LineWidth=2, MarkerColor=8)
plot.Set(h_scale_1, LineColor=9, LineWidth=2, MarkerColor=9)
plot.Set(h_scale_2, LineColor=10, LineWidth=2, MarkerColor=10)
plot.Set(h_scale_3, LineColor=11, LineWidth=2, MarkerColor=11)
plot.Set(h_scale_4, LineColor=12, LineWidth=2, MarkerColor=12)
plot.Set(h_scale_5, LineColor=28, LineWidth=2, MarkerColor=28)

legend.AddEntry(h_env, 'Scale envelope', 'LF')
legend.AddEntry(h_scale_0, '#mu_{R} = 2.0, #mu_{F} = 1.0', 'L')
legend.AddEntry(h_scale_1, '#mu_{R} = 0.5, #mu_{F} = 1.0', 'L')
legend.AddEntry(h_scale_2, '#mu_{R} = 1.0, #mu_{F} = 2.0', 'L')
legend.AddEntry(h_scale_3, '#mu_{R} = 2.0, #mu_{F} = 2.0', 'L')
legend.AddEntry(h_scale_4, '#mu_{R} = 1.0, #mu_{F} = 0.5', 'L')
legend.AddEntry(h_scale_5, '#mu_{R} = 0.5, #mu_{F} = 0.5', 'L')
# h_matrix_fill = h_matrix.Clone()
# plot.Set(h_matrix, LineColor=4, LineWidth=2, MarkerColor=4)
# plot.Set(h_matrix_fill, LineColor=4, LineWidth=2, MarkerColor=4, MarkerSize=0, FillColorAlpha=(4, 0.3))

h_env.Draw('E2SAME')
# h_matrix_fill.Draw('E2SAME')
h_nominal.Draw('HISTSAME')
h_scale_0.Draw('HISTSAME')
h_scale_1.Draw('HISTSAME')
h_scale_2.Draw('HISTSAME')
h_scale_3.Draw('HISTSAME')
h_scale_4.Draw('HISTSAME')
h_scale_5.Draw('HISTSAME')
# h_matrix.Draw('HISTSAME')
h_axes[0].SetMinimum(1E-6)
h_axes[0].SetMaximum(1E0)

# plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)
legend.Draw()
# # # plot.FixBoxPadding(pads[0], legend, 0.05)

# # # # Do the ratio plot
pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].Draw()

# r_MC = plot.MakeRatioHist(h_MC, h_matrix, True, False)
r_env = plot.MakeRatioHist(h_env, h_env, True, False)
r_scale_0 = plot.MakeRatioHist(h_scale_0, h_nominal, True, False)
r_scale_1 = plot.MakeRatioHist(h_scale_1, h_nominal, True, False)
r_scale_2 = plot.MakeRatioHist(h_scale_2, h_nominal, True, False)
r_scale_3 = plot.MakeRatioHist(h_scale_3, h_nominal, True, False)
r_scale_4 = plot.MakeRatioHist(h_scale_4, h_nominal, True, False)
r_scale_5 = plot.MakeRatioHist(h_scale_5, h_nominal, True, False)
r_env.Draw('E2SAME')
# r_matrix_fill.Draw('E2SAME')
r_scale_0.Draw('HISTSAME')
r_scale_1.Draw('HISTSAME')
r_scale_2.Draw('HISTSAME')
r_scale_3.Draw('HISTSAME')
r_scale_4.Draw('HISTSAME')
r_scale_5.Draw('HISTSAME')
plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'Ratio', True, 0.81, 1.19)


# # # Go back and tidy up the axes and frame
pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()


canv.Print('.png')
canv.Print('.pdf')
fout = ROOT.TFile('output_scale_uncertainty.root', 'RECREATE')

NodeToTDir(fout, hists)

fout.Close()
