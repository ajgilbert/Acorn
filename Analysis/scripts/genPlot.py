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

tname = 'WGAnalysis'

samples = {
    'WGSM': 'output/100718/wgamma_2016_v2/wg_gen/WpAToLNuA0j_5f_LO_MLM_pTA_300_0.root',
    'WGTH': 'output/100718/wgamma_2016_v2/wg_gen/Theorists_0.root',
    'WG': 'output/100718/wgamma_2016_v2/wg_gen/WGToMuNuG-EFT-madgraphMLM-stitched_0.root'
}

remap = {
    'WG': 'WGToLNuG-EFT_pTA_300_inf-madgraphMLM'
}


parser = argparse.ArgumentParser()
parser.add_argument('--draw', default='phi1')
parser.add_argument('--abs', action='store_true')
parser.add_argument('--unit-norm', action='store_true')
parser.add_argument('--charge', default='+1')
parser.add_argument('--g_pt', default='300.')
parser.add_argument('--g_pt_max', default='9999999999.')
parser.add_argument('--g_eta', default='3.')
parser.add_argument('--l_pt', default='80.')
parser.add_argument('--l_eta', default='2.4')
parser.add_argument('--n_pt', default='80.')
parser.add_argument('--n_eta', default='9999.')
parser.add_argument('--nparts_max', default='10')
parser.add_argument('--dr', default='3.0')
parser.add_argument('--output', '-o', default='gen_plot')
parser.add_argument('--save-scalings', type=int, default=0, help="1: save absolute 2: save relative")
# parser.add_argument('--ratio', '-o', default='gen_plot')
args = parser.parse_args()

pm_label = 'p' if args.charge == '+1' else 'n' if args.charge == '-1' else 'pn'


fout = ROOT.TFile('%s_w_%s.root' % (args.output, pm_label), 'RECREATE')
hists = Node()

# drawabs = True
# drawvar = 'lhe_phi1'
if args.abs:
    drawvar = 'fabs(%s)' % args.draw
    binning = (5, 0, 3.15)
else:
    drawvar = '%s' % args.draw
    binning = (10, -3.15, 3.15)

pt_bins = [150, 210, 300, 420, 1200]

for name, sa, wt in [
        ('nominal', 'WG', 'wt_C3w_0p0*wt_def'),
        ('C3w_0p1', 'WG', 'wt_C3w_0p1*wt_def'),
        ('C3w_0p2', 'WG', 'wt_C3w_0p2*wt_def'),
        ('C3w_0p4', 'WG', 'wt_C3w_0p4*wt_def'),
        ('C3w_1p0', 'WG', 'wt_C3w_1p0*wt_def'),
        # ('SM', 'WGSM', '1.0'),
        # ('TH', 'WGTH', '1.0')
]:
    if args.charge == '0':
        charge_sel = '1'
    else:
        charge_sel = 'l_charge==%s' % args.charge
    sel = '%s && nparts>=3 && nparts <=%s && g_pt>%s && g_pt<%s && l_pt>%s && n_pt>%s && fabs(l_eta) < %s && fabs(n_eta) < %s && fabs(g_eta) < %s && l_g_dr > %s' % (
        charge_sel, args.nparts_max, args.g_pt, args.g_pt_max, args.l_pt, args.n_pt, args.l_eta, args.n_eta, args.g_eta, args.dr)
    print sel
    hists[name] = Hist('TH1D', binning, sa, [drawvar],
                       sel=sel, wt=wt)
    hists[name + '_2D'] = Hist('TH2D', binning + (len(pt_bins) - 1, array('d', pt_bins)), sa, [drawvar, 'g_pt'], sel=sel, wt=wt)

MultiDraw(hists, samples, tname)


save_scalings = args.save_scalings


if save_scalings >= 1:
    if save_scalings == 2:
        for hname in ['nominal_2D', 'C3w_0p1_2D', 'C3w_0p2_2D', 'C3w_0p4_2D', 'C3w_1p0_2D']:
            htmp = hists[hname]
            for jb in xrange(1, htmp.GetNbinsY() + 1):
                tot = 0.
                for ib in xrange(1, htmp.GetNbinsX() + 1):
                    tot += htmp.GetBinContent(ib, jb)
                for ib in xrange(1, htmp.GetNbinsX() + 1):
                    htmp.SetBinContent(ib, jb, htmp.GetBinContent(ib, jb) / tot)
                    htmp.SetBinError(ib, jb, htmp.GetBinError(ib, jb) / tot)

    for jb in xrange(1, hists['nominal_2D'].GetNbinsY() + 1):
        ymin = hists['nominal_2D'].GetYaxis().GetBinLowEdge(jb)
        ymax = hists['nominal_2D'].GetYaxis().GetBinUpEdge(jb)

        for ib in xrange(1, hists['nominal_2D'].GetNbinsX() + 1):
            xmin = hists['nominal_2D'].GetXaxis().GetBinLowEdge(ib)
            xmax = hists['nominal_2D'].GetXaxis().GetBinUpEdge(ib)
            npoints = 5
            gr = ROOT.TGraphErrors(npoints)
            nom =  hists['nominal_2D'].GetBinContent(ib, jb)
            if nom != 0.:
                gr.SetPoint(0, 0.0, hists['nominal_2D'].GetBinContent(ib, jb) / nom)
                gr.SetPoint(1, 0.1, hists['C3w_0p1_2D'].GetBinContent(ib, jb) / nom)
                gr.SetPoint(2, 0.2, hists['C3w_0p2_2D'].GetBinContent(ib, jb) / nom)
                gr.SetPoint(3, 0.4, hists['C3w_0p4_2D'].GetBinContent(ib, jb) / nom)
                gr.SetPoint(4, 1.0, hists['C3w_1p0_2D'].GetBinContent(ib, jb) / nom)
                gr.SetPointError(0, 0., hists['nominal_2D'].GetBinError(ib, jb) / nom)
                gr.SetPointError(1, 0., hists['C3w_0p1_2D'].GetBinError(ib, jb) / nom)
                gr.SetPointError(2, 0., hists['C3w_0p2_2D'].GetBinError(ib, jb) / nom)
                gr.SetPointError(3, 0., hists['C3w_0p4_2D'].GetBinError(ib, jb) / nom)
                gr.SetPointError(4, 0., hists['C3w_1p0_2D'].GetBinError(ib, jb) / nom)
            canv = ROOT.TCanvas('%s_w_%s_gen_bin_%i_%i' % (args.output, pm_label, jb - 1, ib - 1), '%s_w_%s_gen_bin_%i_%i' % (args.output, pm_label, jb - 1, ib - 1))
            pads = plot.OnePad()
            gr.Draw('APC')
            gr.Print()
            WriteToTFile(gr, fout, '', 'w_%s_gen_bin_%i_%i' % (pm_label, jb - 1, ib - 1))
            plot.DrawTitle(pads[0], '[%g,%g],[%g,%g]' % (ymin, ymax, xmin, xmax), 3)
            canv.Print('.png')
            canv.Print('.pdf')

if args.unit_norm:
    hists.ForEach(lambda x: NormaliseTo(x, 1.0))
# hists.ForEach(lambda x: WidthDivide(x))

canv = ROOT.TCanvas('%s_w_%s' % (args.output, pm_label), '%s_w_%s' % (args.output, pm_label))
pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

# Get the data and create axis hist
h_nominal = hists['nominal']
# h_sm = hists['SM']
# h_th = hists['TH']
h_0p1 = hists['C3w_0p1']
h_0p2 = hists['C3w_0p2']
h_0p4 = hists['C3w_0p4']
h_1p0 = hists['C3w_1p0']

h_axes = [h_nominal.Clone() for x in pads]
for h in h_axes:
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
# legend.AddEntry(h_sm, 'SM', 'L')
# legend.AddEntry(h_th, 'Reference', 'L')
legend.AddEntry(h_0p1, 'C_{3W} = 0.1', 'L')
legend.AddEntry(h_0p2, 'C_{3W} = 0.2', 'L')
legend.AddEntry(h_0p4, 'C_{3W} = 0.4', 'L')
legend.AddEntry(h_1p0, 'C_{3W} = 1.0', 'L')

# plot.Set(h_sm, LineColor=2, LineWidth=2, MarkerColor=2)
# plot.Set(h_th, LineColor=4, LineWidth=2, MarkerColor=4)
plot.Set(h_0p1, LineColor=6, LineWidth=1)
plot.Set(h_0p2, LineColor=7, LineWidth=1)
plot.Set(h_0p4, LineColor=9, LineWidth=1)
plot.Set(h_1p0, LineColor=28, LineWidth=1)

h_nominal.Draw('HISTSAMEE')
# h_sm.Draw('HISTSAMEE')
# h_th.Draw('HISTSAMEE')
h_0p1.Draw('HISTSAME')
h_0p2.Draw('HISTSAME')
h_0p4.Draw('HISTSAME')
h_1p0.Draw('HISTSAME')

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.43)
legend.Draw()
# plot.FixBoxPadding(pads[0], legend, 0.05)

# # Do the ratio plot
pads[1].cd()
pads[1].SetGrid(0, 1)
h_axes[1].Draw()

# r_data = plot.MakeRatioHist(h_data, h_tot, True, False)
# r_nominal = plot.MakeRatioHist(h_nominal, h_nominal, True, False)
r_0p1 = plot.MakeRatioHist(h_0p1, h_nominal, True, False)
r_0p2 = plot.MakeRatioHist(h_0p2, h_nominal, True, False)
r_0p4 = plot.MakeRatioHist(h_0p4, h_nominal, True, False)
r_1p0 = plot.MakeRatioHist(h_1p0, h_nominal, True, False)
r_0p1.Draw('HISTSAME')
r_0p2.Draw('HISTSAME')
r_0p4.Draw('HISTSAME')
r_1p0.Draw('HISTSAME')
# r_data.Draw('SAME')
plot.SetupTwoPadSplitAsRatio(
    pads, plot.GetAxisHist(
        pads[0]), plot.GetAxisHist(pads[1]), 'C_{3W} = X / Nominal', True, 0.64, 1.36)


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
plot.DrawTitle(pads[0], 'W^{%s}#gamma^{} LO - MG5_aMC@NLO 2.4.2' %
               ('+' if args.charge == '+1' else '-' if args.charge == '-1' else '#pm'), 1)

pt_l = ROOT.TPaveText(0.23, 0.65, 0.55, 0.9, 'NDCNB')
pt_l.AddText('Selection: before showering (LHE)')
pt_l.AddText('p_{T}^{#gamma} > %s GeV, |#eta^{#gamma}| < %s' %
             (args.g_pt, args.g_eta))
pt_l.AddText('p_{T}^{l} > %s GeV, |#eta^{l}| < %s' % (args.l_pt, args.l_eta))
pt_l.AddText('p_{T}^{miss} > %s GeV' % args.n_pt)
pt_l.AddText('#DeltaR(l, #gamma) > %s' % args.dr)
pt_l.AddText('nJets (ME) == 0')
plot.Set(pt_l, TextAlign=11, TextFont=42, BorderSize=0)

pt_l.Draw()
# ... and we're done
canv.Print('.png')
canv.Print('.pdf')

for path, name, obj in hists.ListObjects():
    print (path, name, obj)
    WriteToTFile(obj, fout, path, name)

fout.Close()
