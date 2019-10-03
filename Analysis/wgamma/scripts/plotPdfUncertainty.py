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


parser = argparse.ArgumentParser()

parser.add_argument('--task', default='baseline', choices=['eft_region', 'baseline', 'photon_fakes'])
parser.add_argument('--indir', default='/home/files/WGamma_')

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

remap = {
    # 'WG-LO': 'WGToLNuG-madgraphMLM-stitched',
    # 'WG-LO-130': 'WGToLNuG-madgraphMLM-PtG-130',
    # 'WG-LO-500': 'WGToLNuG-madgraphMLM-PtG-500',
    'WG': 'WGToLNuG-amcatnloFXFX'
}

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')


hists = Node()
do_cats = []

# X = SelectionManager()

pt_binning = [30,40,50,60,80,100,200,500,1000]

drawvars = [
    ('gen_p0_pt', pt_binning),
]

npdf = 100

baseline_sel = 'gen_l0_pt>30 && abs(gen_l0_eta)<2.5 && gen_met>0 && gen_p0_pt>30 && abs(gen_p0_eta)<2.5 && gen_l0p0_dr>0.7 && lhe_frixione && gen_met>0'
for var, binning in drawvars:
    for sample in samples:
        hists['inclusive'][var][sample + '_nominal'] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel='%s && (gen_pdgid == 11 || gen_pdgid == 13)' % baseline_sel, wt='wt_def')
        for ipdf in xrange(npdf + 1):
            hists['inclusive'][var][sample + '_%i' % ipdf] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel='%s && (gen_pdgid == 11 || gen_pdgid == 13)' % baseline_sel, wt='wt_def * (1. + wt_pdf_%i)' % ipdf)

MultiDraw(hists, samples, tname, mt_cores=4)


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


h_pdf_env = hists['inclusive']['gen_p0_pt']['WG_0'].Clone()
h_pdf_nlo = hists['inclusive']['gen_p0_pt']['WG_nominal'].Clone()
for ib in xrange(1, h_pdf_env.GetNbinsX() + 1):
    nom = h_pdf_env.GetBinContent(ib)
    dev = 0.
    for ipdf in xrange(1, npdf + 1):
        dev += math.pow(nom - hists['inclusive']['gen_p0_pt']['WG_%i' % ipdf].GetBinContent(ib), 2)
    unc = math.sqrt(dev / float(npdf))
    h_pdf_env.SetBinError(ib, unc)

hists['inclusive']['gen_p0_pt']['WG_0_env'] = h_pdf_env

# matrix_p = "/afs/cern.ch/work/a/agilbert/matrix/MATRIX_v1.0.3/run/ppexnea03_MATRIX/result/run_NNLO_iso_met/gnuplot/histograms/pT_gamma__NNLO_QCD.hist"
# matrix_m = "/afs/cern.ch/work/a/agilbert/matrix/MATRIX_v1.0.3/run/ppenexa03_MATRIX/result/run_NNLO_iso_met/gnuplot/histograms/pT_gamma__NNLO_QCD.hist"
# h_matrix_p = ReadTxtHist(matrix_p)
# h_matrix_m = ReadTxtHist(matrix_m)
# h_matrix_p_sc_lo = ReadTxtHist(matrix_p, column=3)
# h_matrix_p_sc_hi = ReadTxtHist(matrix_p, column=5)
# h_matrix_m_sc_lo = ReadTxtHist(matrix_m, column=3)
# h_matrix_m_sc_hi = ReadTxtHist(matrix_m, column=5)
# h_matrix = h_matrix_p.Clone()
# h_matrix_sc_lo = h_matrix_p_sc_lo.Clone()
# h_matrix_sc_hi = h_matrix_p_sc_hi.Clone()
# h_matrix.Add(h_matrix_m)
# h_matrix_sc_lo.Add(h_matrix_m_sc_lo)
# h_matrix_sc_hi.Add(h_matrix_m_sc_hi)

# h_matrix.Scale(1E-3)
# h_matrix = VariableRebin(h_matrix, pt_binning)
# h_matrix_sc_lo.Scale(1E-3)
# h_matrix_sc_lo = VariableRebin(h_matrix_sc_lo, pt_binning)
# h_matrix_sc_hi.Scale(1E-3)
# h_matrix_sc_hi = VariableRebin(h_matrix_sc_hi, pt_binning)

# print h_matrix.Integral(), h_matrix_sc_lo.Integral(), h_matrix_sc_hi.Integral()

# # h_matrix.Print("range")
# for ib in xrange(1, h_matrix.GetNbinsX() + 1):
#     h_matrix.SetBinError(ib, max(abs(h_matrix.GetBinContent(ib) - h_matrix_sc_lo.GetBinContent(ib)), abs(h_matrix.GetBinContent(ib) - h_matrix_sc_hi.GetBinContent(ib))))
# # h_matrix.Print("range")

# h_matrix.Scale(1., 'width')

# h_MC = hists['inclusive']['gen_p0_pt']['WG-NLO']
# h_MC.Print("range")
# print h_MC.Integral('width')
# h_matrix.Print("range")
# print h_matrix.Integral('width')

canv = ROOT.TCanvas('pdf_uncertainty', 'pdf_uncertainty')
pads = plot.TwoPadSplit(0.35, 0.01, 0.01)

h_axes = [h_pdf_env.Clone() for x in pads]
for h in h_axes:
    h.Reset()

h_axes[0].GetYaxis().SetTitle('a.u.')
h_axes[0].Draw()
pads[0].SetLogy()
# pads[0].SetLogx()
# pads[1].SetLogx()
h_axes[0].SetMinimum(1E-9)

# # # A dict to keep track of the hists
legend = ROOT.TLegend(0.67, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

h_pdf_env_fill = h_pdf_env.Clone()
plot.Set(h_pdf_env, LineColor=2, LineWidth=2, MarkerColor=2)
plot.Set(h_pdf_nlo, LineColor=8, LineWidth=2, MarkerColor=8)
plot.Set(h_pdf_env_fill, LineColor=2, LineWidth=2, MarkerColor=2, MarkerSize=0, FillColorAlpha=(2, 0.3))

legend.AddEntry(h_pdf_env_fill, 'NNPDF31 NNLO', 'LF')
legend.AddEntry(h_pdf_nlo, 'NNPDF30 NLO', 'L')

# h_matrix_fill = h_matrix.Clone()
# plot.Set(h_matrix, LineColor=4, LineWidth=2, MarkerColor=4)
# plot.Set(h_matrix_fill, LineColor=4, LineWidth=2, MarkerColor=4, MarkerSize=0, FillColorAlpha=(4, 0.3))

h_pdf_env_fill.Draw('E2SAME')
# h_matrix_fill.Draw('E2SAME')
h_pdf_env.Draw('HISTSAME')
h_pdf_nlo.Draw('HISTSAME')
# h_matrix.Draw('HISTSAME')

# plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)
legend.Draw()
# # # plot.FixBoxPadding(pads[0], legend, 0.05)

# # # # Do the ratio plot
pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].Draw()

# r_MC = plot.MakeRatioHist(h_MC, h_matrix, True, False)
r_pdf_env_fill = plot.MakeRatioHist(h_pdf_env_fill, h_pdf_env_fill, True, False)
r_pdf_nlo = plot.MakeRatioHist(h_pdf_nlo, h_pdf_env_fill, True, False)
r_pdf_env_fill.Draw('E2SAME')
# r_matrix_fill.Draw('E2SAME')
r_pdf_nlo.Draw('HISTSAME')
plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'ratio to NNLO', True, 0.92, 1.08)


# # # Go back and tidy up the axes and frame
pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()


canv.Print('.png')
canv.Print('.pdf')
fout = ROOT.TFile('output_pdf_uncertainty.root', 'RECREATE')

NodeToTDir(fout, hists)

fout.Close()
