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
    # 'WG': 'output/100718/wgamma_2016_v2/wg_gen/WGToMuNuG-EFT-madgraphMLM-stitched_0.root'
    'WG': '/home/files/190411-gen/wgamma_2016_v2/wg_gen_WGToMuNuG-EFT-madgraphMLM-stitched.root'
}
# /home/files/190411-gen/wgamma_2016_v2/wg_gen_WGToMuNuG-EFT-madgraphMLM-stitched.root
remap = {
    'WG': 'WGToLNuG-EFT_pTA_300_inf-madgraphMLM'
}


# g_pt[]

def GetBinning(binning_str):
    if binning_str.startswith('['):
        binning = [float(v) for v in binning_str[1:-1].split(',')]
        binning = (len(binning) - 1, array('d', binning))
    else:
        binning = binning_str[1:-1].split(',')
        binning = (int(binning[0]), float(binning[1]), float(binning[2]))
    return binning


def ParametrizeBin(x_vals, y_vals, y_val_errs, label, makePlots=True, dropBSM=False, dropInt=False, wsp=None):
    if y_vals[0] == 0.:
        print '>> Skipping bin %s due to zero content' % label
        return
    y_vals_rel = [y / y_vals[0] for y in y_vals]
    y_vals_errs_rel = [y / y_vals[0] for y in y_val_errs]
    npoints = len(x_vals)
    gr = ROOT.TGraphErrors(npoints, array('d', x_vals), array('d', y_vals_rel), array('d', [0] * npoints), array('d', y_vals_errs_rel))
    fn_full = None
    fn_lin = None
    fn_BSM = None
    sig_SM = y_vals_rel[0]
    y_0 = y_vals_rel[0]
    y_0p1 = y_vals_rel[1]
    y_0p2 = y_vals_rel[2]
    sig_BSM = (y_0p2 - 2. * y_0p1 + y_0) / 0.02
    sig_int = (y_0p1 - y_0 - 0.01 * sig_BSM) / 0.1
    print '%.2g, %.2g, %.2g' % (sig_SM, sig_int, sig_BSM)
    if makePlots:
        fn_full = ROOT.TF1("fn_full", "([0] + x*[1] + x*x*[2])/[0]", 0, 1)
        fn_lin = ROOT.TF1("fn_lin", "([0] + x*[1])/[0]", 0, 1)
        fn_BSM = ROOT.TF1("fn_BSM", "([0] + x*x*[1])/[0]", 0, 1)
        fn_full.SetParameter(0, sig_SM)
        fn_full.SetParameter(1, sig_int)
        fn_full.SetParameter(2, sig_BSM)
        fn_lin.SetParameter(0, sig_SM)
        fn_lin.SetParameter(1, sig_int)
        fn_BSM.SetParameter(0, sig_SM)
        fn_BSM.SetParameter(1, sig_BSM)
        canv = ROOT.TCanvas('%s' % (label), '%s' % (label))
        pads = plot.OnePad()
        gr.Draw('APC')
        gr.Print()
        if fn_full is not None:
            fn_full.Draw("SAME")
            fn_lin.SetLineColor(4)
            fn_lin.Draw('SAME')
            fn_BSM.SetLineColor(8)
            fn_BSM.Draw('SAME')
            axis = plot.GetAxisHist(pads[0])
            axis.SetMinimum(0)
            axis.SetMaximum(2)
            axis.GetXaxis().SetRangeUser(0, 1.0)
        # WriteToTFile(gr, fout, '', 'w_%s_gen_bin_%i_%i' % (pm_label, jb - 1, ib - 1))
        # plot.DrawTitle(pads[0], '[%g,%g],[%g,%g]' % (ymin, ymax, xmin, xmax), 3)
        canv.Print('.png')
        canv.Print('.pdf')
    if dropInt:
        sig_int = 0.
    if dropBSM:
        sig_BSM = 0.
    if wsp is not None:
        wsp.factory('expr::%s("%.3g+@0*%.3g+@0*@0*%.3g",c3w[0,0,10])' % (label, sig_SM, sig_int, sig_BSM))
    return (sig_SM, sig_int, sig_BSM)


parser = argparse.ArgumentParser()
parser.add_argument('--draw-x', default='gen_phi', nargs=3)
parser.add_argument('--draw-y', default=None, nargs=3)
parser.add_argument('--unit-norm', action='store_true')
parser.add_argument('--charge', default='+1')
parser.add_argument('--logy', action='store_true')
parser.add_argument('--g_pt', default='300.')
parser.add_argument('--g_pt_max', default='9999999999.')
parser.add_argument('--g_eta', default='3.')
parser.add_argument('--l_pt', default='80.')
parser.add_argument('--l_eta', default='2.4')
parser.add_argument('--n_pt', default='80.')
parser.add_argument('--n_pt_max', default='9999999999.')
parser.add_argument('--n_eta', default='9999.')
parser.add_argument('--nparts_max', default='10')
parser.add_argument('--dr', default='3.0')
parser.add_argument('--output', '-o', default='gen_plot')
parser.add_argument('--ratio-max', default=1.46, type=float)
parser.add_argument('--save-scalings', type=int, default=0, help="1: save absolute 2: save relative")
# parser.add_argument('--ratio', '-o', default='gen_plot')
args = parser.parse_args()

drawvar_x = args.draw_x[0]
binning_x = GetBinning(args.draw_x[1])
label_x = args.draw_x[2]

is2D = False
if args.draw_y is not None:
    is2D = True
    drawvar_y = args.draw_y[0]
    binning_y = GetBinning(args.draw_y[1])


pm_label = 'p' if args.charge == '+1' else 'n' if args.charge == '-1' else 'pn'


fout = ROOT.TFile('%s_w_%s.root' % (args.output, pm_label), 'RECREATE')
wsp = ROOT.RooWorkspace('w','w')
hists = Node()



# drawabs = True
# drawvar = 'lhe_phi1'
# pt_bins = [150, 210, 300, 420, 600, 850, 1200]
# if args.abs:
#     drawvar = 'fabs(%s)' % args.draw
#     binning = (5, 0, 3.15)
# else:
#     drawvar = '%s' % args.draw
#     if 'phi' in drawvar:
#         binning = (10, -3.15, 3.15)
#     elif 'g_pt' in drawvar:
#         binning = (20, 0, 1000)

compute = [
    ('X_0p1', 0.1),
    ('X_0p2', 0.2),
    ('X_0p4', 0.4),
    ('X_0p67', 0.67),
    ('X_1p0', 1.0),
]

for name, sa, wt in [
        ('nominal', 'WG', 'wt_C3w_0p0*wt_def'),
        ('C3w_0p1', 'WG', 'wt_C3w_0p1*wt_def'),
        ('C3w_0p2', 'WG', 'wt_C3w_0p2*wt_def'),
        ('C3w_0p4', 'WG', 'wt_C3w_0p4*wt_def'),
        ('C3w_0p67', 'WG', 'wt_C3w_0p67*wt_def'),
        ('C3w_1p0', 'WG', 'wt_C3w_1p0*wt_def'),
        # ('SM', 'WGSM', '1.0'),
        # ('TH', 'WGTH', '1.0')
]:
    if args.charge == '0':
        charge_sel = '1'
    else:
        charge_sel = 'l_charge==%s' % args.charge
    sel = '%s && nparts>=3 && nparts <=%s && g_pt>%s && g_pt<%s && l_pt>%s && n_pt>%s && n_pt<%s && fabs(l_eta) < %s && fabs(n_eta) < %s && fabs(g_eta) < %s && l_g_dr > %s' % (
        charge_sel, args.nparts_max, args.g_pt, args.g_pt_max, args.l_pt, args.n_pt, args.n_pt_max, args.l_eta, args.n_eta, args.g_eta, args.dr)
    print sel

    if is2D:
        hists[name] = Hist('TH2D', binning_x + binning_y, sa, [drawvar_x, drawvar_y], sel=sel, wt=wt)
    else:
        hists[name] = Hist('TH1D', binning_x, sa, [drawvar_x], sel=sel, wt=wt)
    hists[name + '_phi_reco'] = Hist('TH2D', (40, -3.15, 3.15, 40, -3.15, 3.15), sa, ['true_phi', 'gen_phi'], sel=sel, wt=wt)

MultiDraw(hists, samples, tname)

for label, val in compute:
    hists[label] = hists['nominal'].Clone()

x_vals = [0., 0.1, 0.2, 0.4, 0.67, 1.0]
y_labels= ['nominal', 'C3w_0p1', 'C3w_0p2', 'C3w_0p4', 'C3w_0p67', 'C3w_1p0']
bin_scalings = []
if not is2D:
    for ib in xrange(1, hists['nominal'].GetNbinsX() + 1):
        y_vals = [hists[h].GetBinContent(ib) for h in y_labels]
        y_vals_err = [hists[h].GetBinError(ib) for h in y_labels]
        bin_scalings.append(ParametrizeBin(x_vals, y_vals, y_vals_err, '%s_w_%s_bin_%i' % (args.output, pm_label, ib), wsp=wsp))
        scale = bin_scalings[ib - 1]
        for label, val in compute:
            scale_factor = (scale[0] + val * scale[1] + val * val * scale[2])
            hists[label].SetBinContent(ib, hists[label].GetBinContent(ib) * scale_factor)

else:
    for ib in xrange(1, hists['nominal'].GetNbinsX() + 1):
        bin_scalings.append(list())
        for jb in xrange(1, hists['nominal'].GetNbinsY() + 1):
            y_vals = [hists[h].GetBinContent(ib, jb) for h in y_labels]
            y_vals_err = [hists[h].GetBinError(ib, jb) for h in y_labels]
            bin_scalings[ib - 1].append(ParametrizeBin(x_vals, y_vals, y_vals_err, '%s_w_%s_bin_%i_%i' % (args.output, pm_label, ib, jb), wsp=wsp))
            scale = bin_scalings[ib - 1][jb - 1]
            for label, val in compute:
                scale_factor = (scale[0] + val * scale[1] + val * val * scale[2])
                hists[label].SetBinContent(ib, jb, hists[label].GetBinContent(ib, jb) * scale_factor)


if not is2D:
    if args.unit_norm:
        hists.ForEach(lambda x: NormaliseTo(x, 1.0))
    hists.ForEach(lambda x: WidthDivide(x))

    canv = ROOT.TCanvas('%s_w_%s' % (args.output, pm_label), '%s_w_%s' % (args.output, pm_label))
    pads = plot.TwoPadSplit(0.27, 0.01, 0.01)

    # Get the data and create axis hist
    h_nominal = hists['nominal']
    # h_sm = hists['SM']
    # h_th = hists['TH']
    h_0p1 = hists['X_0p1']
    h_0p2 = hists['X_0p2']
    h_0p4 = hists['X_0p4']
    h_1p0 = hists['X_1p0']

    h_axes = [h_nominal.Clone() for x in pads]
    for h in h_axes:
        h.Reset()


    # if 'true_phi' in drawvar:
    #     labelvar = '#varphi_{true}'
    # elif 'g_pt' in drawvar:
    #     labelvar = 'p_{T}^{#gamma} (GeV)'
    # else:
    #     labelvar = '#varphi_{gen}'
    # if args.abs:
    #     labelvar = '|%s|' % labelvar
    h_axes[1].GetXaxis().SetTitle(label_x)

    h_axes[0].GetYaxis().SetTitle('a.u.')
    h_axes[0].Draw()
    if args.logy:
        pads[0].SetLogy()
        h_axes[0].SetMinimum(0.001)

    # A dict to keep track of the hists
    legend = ROOT.TLegend(0.67, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

    legend.AddEntry(h_nominal, 'C_{3W} = 0.0 (EWDim6)', 'L')
    # legend.AddEntry(h_sm, 'SM', 'L')
    # legend.AddEntry(h_th, 'Reference', 'L')
    legend.AddEntry(h_0p1, 'C_{3W} = 0.1', 'L')
    legend.AddEntry(h_0p2, 'C_{3W} = 0.2', 'L')
    legend.AddEntry(h_0p4, 'C_{3W} = 0.4', 'L')
    legend.AddEntry(h_1p0, 'C_{3W} = 1.0', 'L')

    # plot.Set(h_sm, LineColor=2, LineWidth=2, MarkerColor=2)
    # plot.Set(h_th, LineColor=4, LineWidth=2, MarkerColor=4)
    plot.Set(h_nominal, LineWidth=2)
    plot.Set(h_0p1, LineColor=6, LineWidth=2)
    plot.Set(h_0p2, LineColor=ROOT.kGreen-3, LineWidth=2)
    plot.Set(h_0p4, LineColor=9, LineWidth=2)
    plot.Set(h_1p0, LineColor=28, LineWidth=2)

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
            pads[0]), plot.GetAxisHist(pads[1]), 'C_{3W} = X / Nominal', True, 0.54, args.ratio_max)


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
    pt_l.AddText('Selection:')
    pt_l.AddText('p_{T}^{#gamma} > %s GeV, |#eta^{#gamma}| < %s' %
                 (args.g_pt, args.g_eta))
    pt_l.AddText('p_{T}^{l} > %s GeV, |#eta^{l}| < %s' % (args.l_pt, args.l_eta))
    pt_l.AddText('p_{T}^{miss} > %s GeV' % args.n_pt)
    pt_l.AddText('#DeltaR(l, #gamma) > %s' % args.dr)
    # pt_l.AddText('nJets (ME) >= 0')
    plot.Set(pt_l, TextAlign=11, TextFont=42, BorderSize=0)

    pt_l.Draw()
    # ... and we're done
    canv.Print('.png')
    canv.Print('.pdf')


# save_scalings = args.save_scalings


# if save_scalings >= 1:
#     if save_scalings == 2:
#         for hname in ['nominal_2D', 'C3w_0p1_2D', 'C3w_0p2_2D', 'C3w_0p4_2D', 'C3w_1p0_2D']:
#             htmp = hists[hname]
#             for jb in xrange(1, htmp.GetNbinsY() + 1):
#                 tot = 0.
#                 for ib in xrange(1, htmp.GetNbinsX() + 1):
#                     tot += htmp.GetBinContent(ib, jb)
#                 for ib in xrange(1, htmp.GetNbinsX() + 1):
#                     htmp.SetBinContent(ib, jb, htmp.GetBinContent(ib, jb) / tot)
#                     htmp.SetBinError(ib, jb, htmp.GetBinError(ib, jb) / tot)

#     for jb in xrange(1, hists['nominal_2D'].GetNbinsY() + 1):
#         ymin = hists['nominal_2D'].GetYaxis().GetBinLowEdge(jb)
#         ymax = hists['nominal_2D'].GetYaxis().GetBinUpEdge(jb)

#         for ib in xrange(1, hists['nominal_2D'].GetNbinsX() + 1):
#             xmin = hists['nominal_2D'].GetXaxis().GetBinLowEdge(ib)
#             xmax = hists['nominal_2D'].GetXaxis().GetBinUpEdge(ib)
#             npoints = 6
#             gr = ROOT.TGraphErrors(npoints)
#             nom = hists['nominal_2D'].GetBinContent(ib, jb)
#             if nom != 0.:
#                 gr.SetPoint(0, 0.0, hists['nominal_2D'].GetBinContent(ib, jb) / nom)
#                 gr.SetPoint(1, 0.1, hists['C3w_0p1_2D'].GetBinContent(ib, jb) / nom)
#                 gr.SetPoint(2, 0.2, hists['C3w_0p2_2D'].GetBinContent(ib, jb) / nom)
#                 gr.SetPoint(3, 0.4, hists['C3w_0p4_2D'].GetBinContent(ib, jb) / nom)
#                 gr.SetPoint(4, 0.67, hists['C3w_0p67_2D'].GetBinContent(ib, jb) / nom)
#                 gr.SetPoint(5, 1.0, hists['C3w_1p0_2D'].GetBinContent(ib, jb) / nom)
#                 gr.SetPointError(0, 0., hists['nominal_2D'].GetBinError(ib, jb) / nom)
#                 gr.SetPointError(1, 0., hists['C3w_0p1_2D'].GetBinError(ib, jb) / nom)
#                 gr.SetPointError(2, 0., hists['C3w_0p2_2D'].GetBinError(ib, jb) / nom)
#                 gr.SetPointError(3, 0., hists['C3w_0p4_2D'].GetBinError(ib, jb) / nom)
#                 gr.SetPointError(4, 0., hists['C3w_0p67_2D'].GetBinError(ib, jb) / nom)
#                 gr.SetPointError(5, 0., hists['C3w_1p0_2D'].GetBinError(ib, jb) / nom)
#             canv = ROOT.TCanvas('%s_w_%s_gen_bin_%i_%i' % (args.output, pm_label, jb - 1, ib - 1), '%s_w_%s_gen_bin_%i_%i' % (args.output, pm_label, jb - 1, ib - 1))
#             pads = plot.OnePad()
#             gr.Draw('APC')
#             gr.Print()
#             WriteToTFile(gr, fout, '', 'w_%s_gen_bin_%i_%i' % (pm_label, jb - 1, ib - 1))
#             plot.DrawTitle(pads[0], '[%g,%g],[%g,%g]' % (ymin, ymax, xmin, xmax), 3)
#             canv.Print('.png')
#             canv.Print('.pdf')



for path, name, obj in hists.ListObjects():
    print (path, name, obj)
    WriteToTFile(obj, fout, path, name)
wsp.Write()
fout.Close()
