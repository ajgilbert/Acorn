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
parser.add_argument('--indir', default='root://eoscms.cern.ch//store/cmst3/user/agilbert/190630-full/wgamma_2016_v4/WGamma_')

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

remap = {
    # 'WG-LO': 'WGToLNuG-madgraphMLM-stitched',
    # 'WG-LO-130': 'WGToLNuG-madgraphMLM-PtG-130',
    # 'WG-LO-500': 'WGToLNuG-madgraphMLM-PtG-500',
    'WG-NLO': 'WGToLNuG-amcatnloFXFX-stitched'
}

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')


hists = Node()
do_cats = []

# X = SelectionManager()

pt_binning = [0, 15, 30, 45, 60, 90, 130, 200, 300, 400, 500, 1000]

drawvars = [
    ('lhe_p0_pt', pt_binning),
    ('gen_p0_pt', pt_binning),
    ('lhe_l0_pt', pt_binning),
    ('gen_l0_pt', pt_binning),
    ('lhe_p0_eta', (50, -5, 5)),
    ('gen_p0_eta', (50, -5, 5)),
    ('lhe_l0_eta', (50, -5, 5)),
    ('gen_l0_eta', (50, -5, 5)),
]

# This is the selection used by the MATRIX authors:
# baseline_sel = 'gen_p0_pt > 15 && abs(gen_p0_eta) < 2.5 && gen_l0_pt > 25 && abs(gen_l0_eta) < 2.5 && gen_met > 35 && gen_m0p0_dr > 0.7'
# baseline_sel = 'gen_p0_pt > 15 && abs(gen_p0_eta) < 2.6 && gen_l0_pt > 15 && abs(gen_l0_eta) < 999.0'
baseline_sel = 'gen_l0_q==+1 && gen_l0_pt>35 && abs(gen_l0_eta)<2.4 && gen_met>30 && gen_p0_pt>30 && (gen_p0_eta)<2.5 && gen_l0p0_dr>0.7 && lhe_frixione'
for var, binning in drawvars:
    for sample in samples:
        hists['inclusive'][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel='%s && gen_pdgid == 11' % baseline_sel, wt='wt_def')

MultiDraw(hists, samples, tname, mt_cores=4)


xsecs = {
    'WG-LO': 465.4,
    'WG-NLO': 178.9
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

matrix_p = "/afs/cern.ch/work/a/agilbert/matrix/MATRIX_v1.0.3/run/ppexnea03_MATRIX/result/run_test_NNLO_iso_met/gnuplot/histograms/pT_gamma__NNLO_QCD.hist"
h_matrix_NLO_p = ReadTxtHist(matrix_p)
h_matrix_NLO_m = ReadTxtHist(matrix_p)
h_matrix_p_sc_lo = ReadTxtHist(matrix_p, column=3)
h_matrix_p_sc_hi = ReadTxtHist(matrix_p, column=5)
h_matrix_NLO = h_matrix_NLO_p.Clone()
# h_matrix_NLO.Add(h_matrix_NLO_m)

h_matrix_NLO.Scale(1E-3)
h_matrix_NLO = VariableRebin(h_matrix_NLO, pt_binning)
h_matrix_p_sc_lo.Scale(1E-3)
h_matrix_p_sc_lo = VariableRebin(h_matrix_p_sc_lo, pt_binning)
h_matrix_p_sc_hi.Scale(1E-3)
h_matrix_p_sc_hi = VariableRebin(h_matrix_p_sc_hi, pt_binning)
for ib in xrange(1, h_matrix_NLO.GetNbinsX() + 1):
    h_matrix_NLO.SetBinError(ib, max(abs(h_matrix_NLO.GetBinContent(ib) - h_matrix_p_sc_lo.GetBinContent(ib)), abs(h_matrix_NLO.GetBinContent(ib) - h_matrix_p_sc_hi.GetBinContent(ib))))
h_matrix_NLO.Scale(1., 'width')

h_MC = hists['inclusive']['gen_p0_pt']['WG-NLO']
h_MC.Print("range")
print h_MC.Integral('width')
h_matrix_NLO.Print("range")
print h_matrix_NLO.Integral('width')

#h_matrix_NLO.Divide(h_MC)
#h_matrix_NLO.Print("range")

canv = ROOT.TCanvas('test', 'test')
pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

# # Get the data and create axis hist
# var = 'gen_p0_pt'
# h_LO_inc = hists['inclusive'][var]['WG-LO-inc']
# # h_LO_130 = hists['inclusive'][var]['WG-LO-130']
# h_LO_inc_cut_130 = hists['inclusive'][var]['WG-LO-inc-cut-130']
# # h_LO_500 = hists['inclusive'][var]['WG-LO-500']
# h_NLO_inc = hists['inclusive'][var]['WG-NLO-inc']

h_axes = [h_MC.Clone() for x in pads]
for h in h_axes:
    h.Reset()

# h_LO_inc.Scale(1. / h_LO_inc.Integral(), 'width')
# h_NLO_inc.Scale(1. / h_NLO_inc.Integral(), 'width')
# # h_LO_130.Scale(frac_130 / h_LO_130.Integral(), 'width')
# h_LO_inc_cut_130.Scale(frac_130 / h_LO_inc_cut_130.Integral(), 'width')
# # h_LO_500.Scale(frac_500 / h_LO_500.Integral(), 'width')

# # if 'true_phi' in drawvar:
# #     labelvar = '#varphi_{true}'
# # elif 'g_pt' in drawvar:
# #     labelvar = 'p_{T}^{#gamma} (GeV)'
# # else:
# #     labelvar = '#varphi_{gen}'
# # if args.abs:
# #     labelvar = '|%s|' % labelvar
# # h_axes[1].GetXaxis().SetTitle(labelvar)

h_axes[0].GetYaxis().SetTitle('a.u.')
h_axes[0].Draw()
pads[0].SetLogy()
h_axes[0].SetMinimum(1E-7)

# # A dict to keep track of the hists
legend = ROOT.TLegend(0.67, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

legend.AddEntry(h_MC, 'MG NLO inclusive', 'L')
# # legend.AddEntry(h_sm, 'SM', 'L')
# # legend.AddEntry(h_th, 'Reference', 'L')
legend.AddEntry(h_matrix_NLO, 'MATRIX NNLO', 'L')
# # legend.AddEntry(h_LO_130, 'LO 130', 'L')
# # legend.AddEntry(h_LO_inc_cut_130, 'LO inc., 130 cut', 'L')
# # legend.AddEntry(h_LO_500, 'LO 500', 'L')

plot.Set(h_MC, LineColor=2, LineWidth=2, MarkerColor=2)
plot.Set(h_matrix_NLO, LineColor=4, LineWidth=2, MarkerColor=4)
# plot.Set(h_LO_inc, LineColor=1, LineWidth=2, MarkerColor=1, MarkerSize=0.5)
# plot.Set(h_NLO_inc, LineColor=ROOT.kGreen-3, LineWidth=1, MarkerColor=ROOT.kGreen-3, MarkerSize=0.5)
# # plot.Set(h_LO_130, LineColor=9, LineWidth=1)
# # plot.Set(h_LO_inc_cut_130, LineColor=2, LineWidth=2)
# # plot.Set(h_LO_500, LineColor=28, LineWidth=1)

h_MC.Draw('HISTSAMEE')
h_matrix_NLO.Draw('HISTSAMEE2')
# # h_th.Draw('HISTSAMEE')
# h_NLO_inc.Draw('HISTSAMEE')
# # h_LO_130.Draw('HISTSAME')
# # h_LO_inc_cut_130.Draw('HISTSAME')
# # h_LO_500.Draw('HISTSAME')

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)
legend.Draw()
# # plot.FixBoxPadding(pads[0], legend, 0.05)

# # # Do the ratio plot
pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].Draw()

r_MC = plot.MakeRatioHist(h_MC, h_matrix_NLO, True, True)
# # r_nominal = plot.MakeRatioHist(h_nominal, h_nominal, True, False)
# r_NLO_inc = plot.MakeRatioHist(h_NLO_inc, h_LO_inc, True, True)
# # r_LO_130 = plot.MakeRatioHist(h_LO_130, h_LO_inc, True, True)
# # r_LO_inc_cut_130 = plot.MakeRatioHist(h_LO_130, h_LO_inc, True, True)
# # r_LO_500 = plot.MakeRatioHist(h_LO_inc_cut_130, h_LO_inc, True, True)
r_MC.Draw('SAMEE')
# # r_LO_130.Draw('SAMEE')
# # r_LO_inc_cut_130.Draw('SAMEE')
# # r_LO_500.Draw('SAMEE')
# # r_data.Draw('SAME')
plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'MC / NNLO', True, 0.8, 1.2)


# # Go back and tidy up the axes and frame
pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

# # CMS logo
# # plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)
# # plot.DrawTitle(pads[0], '0.1 fb^{-1} (13 TeV)', 3)

# # latex = ROOT.TLatex()
# # plot.Set(latex, NDC=None, TextFont=42, TextSize=0.03)
# # # latex.DrawLatex(0.20, 0.75, args.title)
# # plot.DrawTitle(pads[0], '#sqrt{s} = 13 TeV', 3)
# # plot.DrawTitle(pads[0], 'W^{%s}#gamma^{} LO - MG5_aMC@NLO 2.4.2' %
# #                ('+' if args.charge == '+1' else '-' if args.charge == '-1' else '#pm'), 1)

# # pt_l = ROOT.TPaveText(0.23, 0.65, 0.55, 0.9, 'NDCNB')
# # pt_l.AddText('Selection:')
# # pt_l.AddText('p_{T}^{#gamma} > %s GeV, |#eta^{#gamma}| < %s' %
# #              (args.g_pt, args.g_eta))
# # pt_l.AddText('p_{T}^{l} > %s GeV, |#eta^{l}| < %s' % (args.l_pt, args.l_eta))
# # pt_l.AddText('p_{T}^{miss} > %s GeV' % args.n_pt)
# # pt_l.AddText('#DeltaR(l, #gamma) > %s' % args.dr)
# # # pt_l.AddText('nJets (ME) >= 0')
# # plot.Set(pt_l, TextAlign=11, TextFont=42, BorderSize=0)

# # pt_l.Draw()
# # ... and we're done
canv.Print('.png')
canv.Print('.pdf')




# # with open('input/cfg_wgamma_2016_v2.json') as jsonfile:
# #     cfg = json.load(jsonfile)
# #     sample_cfg = cfg['samples']

# # for sample in samples:
# #     f = ROOT.TFile(samples[sample])
# #     sample_cfg[remap[sample]]['events'] = f.Get('counters').GetBinContent(2)
# #     f.Close()

# # for path, hname, obj in hists.ListObjects():
# #     name = obj.sample
# #     if name is not 'data_obs':
# #         tgt_lumi = sample_cfg[remap['data_obs']]['lumi']
# #         events = sample_cfg[remap[name]]['events']
# #         xsec = sample_cfg[remap[name]]['xsec']
# #         scale = tgt_lumi * xsec / events
# #         obj.Scale(scale)

# fout = ROOT.TFile('output_wgamma_signal_%s.root' % args.task, 'RECREATE')

# # for path, node in hists.ListNodes(withObjects=True):
# #     node['VV'] = sum([node[X] for X in dibosons[1:]], node[dibosons[0]])
# #     node['VV_F'] = sum([node['%s_F' % X] for X in dibosons[1:]], node['%s_F' % dibosons[0]])
# #     node['VV_R'] = sum([node['%s_R' % X] for X in dibosons[1:]], node['%s_R' % dibosons[0]])
# #     node['Total_R'] = node['WG'] + node['TT_R'] + node['DY_R'] + node['VV_R']
# #     node['Total_F'] = node['W_F'] + node['TT_F'] + node['DY_F'] + node['VV_F']

# NodeToTDir(fout, hists)

# fout.Close()
