import ROOT
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
import argparse
import math

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

tname = 'WGAnalysis'


def ParametrizeBin(x_vals, y_vals, y_val_errs, sumw2, sA, sB, sA_2, sB_2, label, makePlots=False, dropBSM=False, dropInt=False, wsp=None, binStr=''):
    if y_vals[0] == 0.:
        print '>> Skipping bin %s due to zero content' % label
        return
    y_vals_rel = [y / y_vals[0] for y in y_vals]
    # y_vals_errs_rel = [y / y_vals[0] for y in y_val_errs]
    y_vals_errs_rel = [0] * len(y_val_errs)
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
    A = sA / y_vals[0]
    B = sB / y_vals[0]

    neff = math.pow(y_vals[0], 2) / sumw2[0]
    print 'neff: %g' % neff
    stddev2_A = (sA_2 / y_vals[0]) - A * A
    stddev2_B = (sB_2 / y_vals[0]) - B * B
    if stddev2_A < 0.:
        print 'Error, stddev2_A = %g' % stddev2_A
        stddev_A = 0.
    else:
        stddev_A = math.sqrt(stddev2_A)
    if stddev2_B < 0.:
        print 'Error, stddev2_B = %g' % stddev2_B
        stddev_B = 0.
    else:
        stddev_B = math.sqrt(stddev2_B)
    stderr_A = stddev_A / math.sqrt(neff)
    stderr_B = stddev_B / math.sqrt(neff)

    print 'Nom: %.2g, %.2g, %.2g' % (sig_SM, sig_int, sig_BSM)
    print 'Alt: %.2g, %.2g +/- %.2g, %.2g +/- %2.g' % (sig_SM, A, stderr_A, B, stderr_B)
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
        plot.Set(gr, MarkerColor=2, LineColor=2, LineWidth=2)
        pads = plot.OnePad()
        gr.Draw('AP')
        gr.Print()
        if fn_full is not None:
            plot.Set(fn_full, MarkerColor=2, LineColor=2, LineWidth=2)
            fn_full.Draw("SAME")
            plot.Set(fn_lin, MarkerColor=4, LineColor=4, LineWidth=2)
            fn_lin.Draw('SAME')
            plot.Set(fn_BSM, MarkerColor=1, LineColor=1, LineWidth=2)
            fn_BSM.Draw('SAME')
            axis = plot.GetAxisHist(pads[0])
            axis.SetMinimum(min(axis.GetMinimum(), fn_lin.Eval(1)))
            axis.SetMaximum(max(axis.GetMaximum(), fn_BSM.Eval(1)))
            # axis.SetMaximum(2)
            plot.Set(axis.GetXaxis(), Title='C_{3W} (TeV^{-2})')
            plot.Set(axis.GetYaxis(), Title='#sigma(C_{3W})/#sigma_{SM}')
            legend = ROOT.TLegend(0.2, 0.86 - 0.04 * 3, 0.4, 0.91, '', 'NBNDC')
            legend.AddEntry(fn_full, 'Full', 'L')
            legend.AddEntry(fn_lin, 'Int. only', 'L')
            legend.AddEntry(fn_BSM, 'BSM^2 only', 'L')
            legend.Draw()
            # axis.GetXaxis().SetRangeUser(0, 1.0)
        # WriteToTFile(gr, fout, '', 'w_%s_gen_bin_%i_%i' % (pm_label, jb - 1, ib - 1))
        plot.DrawTitle(pads[0], '#scale[0.7]{#sigma/#sigma_{SM} = 1 + (%.2f)#upointC_{3W} + (%.2f)#upointC_{3W}^{2}}' % (sig_int, sig_BSM), 1)
        plot.DrawTitle(pads[0], '#scale[0.7]{%s}' % binStr, 3)
        canv.Print(plotdir + '/%s.png' % label)
        canv.Print(plotdir + '/%s.pdf' % label)
    if dropInt:
        sig_int = 0.
    if dropBSM:
        sig_BSM = 0.
    if wsp is not None:
        wsp.factory('expr::%s("%.3g+@0*%.3g+@1*@0*@0*%.3g",c3w[0,0,10],withBSM[1])' % (label, sig_SM, sig_int, sig_BSM))
    return (sig_SM, sig_int, sig_BSM, stderr_A, stderr_B)


parser = argparse.ArgumentParser()
parser.add_argument('--sample', default='')
parser.add_argument('--pre-wt', default='wt')
parser.add_argument('--pre-var', default='gen')
parser.add_argument('--draw-x', default='gen_phi', nargs=3)
parser.add_argument('--draw-y', default=None, nargs=3)
parser.add_argument('--unit-norm', action='store_true')
parser.add_argument('--charge', default='+1')
parser.add_argument('--logy', action='store_true')
parser.add_argument('--int-only', action='store_true')
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
parser.add_argument('--label', default='gen_plot')
parser.add_argument('--ratio-max', default=1.46, type=float)
parser.add_argument('--save-scalings', type=int, default=0, help="1: save absolute 2: save relative")
parser.add_argument('--plot-dir', default='.')
args = parser.parse_args()

plotdir = args.plot_dir

drawvar_x = args.draw_x[0]
binning_x = BinningFromStr(args.draw_x[1])
label_x = args.draw_x[2]

is2D = False
if args.draw_y is not None:
    is2D = True
    drawvar_y = args.draw_y[0]
    binning_y = BinningFromStr(args.draw_y[1])
    label_y = args.draw_y[2]

pm_label = 'p' if args.charge == '+1' else 'n' if args.charge == '-1' else 'x'


fout = ROOT.TFile('%s_%s_%s.root' % (args.output, args.label, pm_label), 'RECREATE')
wsp = ROOT.RooWorkspace('w','w')
hists = Node()

compute = [
    ('X_0p1', 0.1),
    ('X_0p2', 0.2),
    ('X_0p4', 0.4),
    ('X_0p67', 0.67),
    ('X_1p0', 1.0),
]

pre_wt = args.pre_wt
gen = args.pre_var


for name, sa, wt in [
        ('nominal', 'WG', '%s_C3w_0p0 * wt_def' % pre_wt),
        ('A', 'WG', '(20*(%s_C3w_0p1/%s_C3w_0p0) - 5*(%s_C3w_0p2/%s_C3w_0p0) - 15) * %s_C3w_0p0 * wt_def' % (pre_wt, pre_wt, pre_wt, pre_wt, pre_wt)),
        ('B', 'WG', '(-100*(%s_C3w_0p1/%s_C3w_0p0) + 50*(%s_C3w_0p2/%s_C3w_0p0) + 50) * %s_C3w_0p0 * wt_def' % (pre_wt, pre_wt, pre_wt, pre_wt, pre_wt)),
        ('A_2', 'WG', '(20*(%s_C3w_0p1/%s_C3w_0p0) - 5*(%s_C3w_0p2/%s_C3w_0p0) - 15) * (20*(%s_C3w_0p1/%s_C3w_0p0) - 5*(%s_C3w_0p2/%s_C3w_0p0) - 15) * %s_C3w_0p0  * wt_def' % (pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt)),
        ('B_2', 'WG', '(-100*(%s_C3w_0p1/%s_C3w_0p0) + 50*(%s_C3w_0p2/%s_C3w_0p0) + 50) * (-100*(%s_C3w_0p1/%s_C3w_0p0) + 50*(%s_C3w_0p2/%s_C3w_0p0) + 50) * %s_C3w_0p0  * wt_def' % (pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt, pre_wt)),
        ('C3w_0p1', 'WG', '%s_C3w_0p1 * wt_def' % pre_wt),
        ('C3w_0p2', 'WG', '%s_C3w_0p2 * wt_def' % pre_wt),
        ('C3w_0p4', 'WG', '%s_C3w_0p4 * wt_def' % pre_wt),
        ('C3w_0p67', 'WG', '%s_C3w_0p67 * wt_def' % pre_wt),
        ('C3w_1p0', 'WG', '%s_C3w_1p0 * wt_def' % pre_wt),
        ('nominal_2', 'WG', '%s_C3w_0p0 * %s_C3w_0p0 * wt_def' % (pre_wt, pre_wt)),
        ('C3w_0p1_2', 'WG', '%s_C3w_0p1 * %s_C3w_0p1 * wt_def' % (pre_wt, pre_wt)),
        ('C3w_0p2_2', 'WG', '%s_C3w_0p2 * %s_C3w_0p2 * wt_def' % (pre_wt, pre_wt)),
        ('C3w_0p4_2', 'WG', '%s_C3w_0p4 * %s_C3w_0p4 * wt_def' % (pre_wt, pre_wt)),
        ('C3w_0p67_2', 'WG', '%s_C3w_0p67 * %s_C3w_0p67 * wt_def' % (pre_wt, pre_wt)),
        ('C3w_1p0_2', 'WG', '%s_C3w_1p0 * %s_C3w_1p0 * wt_def' % (pre_wt, pre_wt)),
]:
    if args.charge == '0':
        charge_sel = '1'
    else:
        charge_sel = '%s_l0_q==%s' % (gen, args.charge)
    sel = '%s_pdgid==13 && %s && nparts>=1 && nparts <=%s && %s_p0_pt>%s && %s_p0_pt<%s && %s_l0_pt>%s && %s_n0_pt>%s && %s_n0_pt<%s && fabs(%s_l0_eta) < %s && fabs(%s_n0_eta) < %s && fabs(%s_p0_eta) < %s && %s_l0p0_dr > %s && %s_p0_frixione' % (
        gen, charge_sel, args.nparts_max, gen, args.g_pt, gen, args.g_pt_max, gen, args.l_pt, gen, args.n_pt, gen, args.n_pt_max, gen, args.l_eta, gen, args.n_eta, gen, args.g_eta, gen, args.dr, gen)
    print sel

    if is2D:
        hists[name] = Hist('TH2D', binning_x + binning_y, sa, [drawvar_x, drawvar_y], sel=sel, wt=wt)
    else:
        hists[name] = Hist('TH1D', binning_x, sa, [drawvar_x], sel=sel, wt=wt)
    hists[name + '_phi_reco'] = Hist('TH2D', (40, -3.15, 3.15, 40, -3.15, 3.15), sa, ['gen_true_phi', 'gen_phi'], sel=sel, wt=wt)

samples = {
    'WG': args.sample
}
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
        sumw2 = [hists[h].GetSumw2()[ib] for h in y_labels]
        sumwx2 = [hists['%s_2' % h].GetBinContent(ib) for h in y_labels]
        A = hists['A'].GetBinContent(ib)
        B = hists['B'].GetBinContent(ib)
        A_2 = hists['A_2'].GetBinContent(ib)
        B_2 = hists['B_2'].GetBinContent(ib)
        bin_str = '%g #leq %s < %g' % (hists[h].GetXaxis().GetBinLowEdge(ib), label_x, hists[h].GetXaxis().GetBinUpEdge(ib))
        bin_scalings.append(ParametrizeBin(x_vals, y_vals, y_vals_err, sumw2, A, B, A_2, B_2, '%s_%s_%i' % (args.label, pm_label, ib - 1), wsp=wsp, binStr=bin_str))
        scale = bin_scalings[ib - 1]
        for label, val in compute:
            if args.int_only:
                scale_factor = (scale[0] + val * scale[1])
            else:
                scale_factor = (scale[0] + val * scale[1] + val * val * scale[2])
            hists[label].SetBinContent(ib, hists[label].GetBinContent(ib) * scale_factor)

    for ib in xrange(0, hists['nominal'].GetNbinsX()):
        print '%-10i %-20.2f %15.2f +/- %-15.2f  %15.2f +/- %-15.2f' % (ib, bin_scalings[ib][0], bin_scalings[ib][1], bin_scalings[ib][3], bin_scalings[ib][2], bin_scalings[ib][4])

else:
    for ib in xrange(1, hists['nominal'].GetNbinsX() + 1):
        bin_scalings.append(list())
        for jb in xrange(1, hists['nominal'].GetNbinsY() + 1):
            y_vals = [hists[h].GetBinContent(ib, jb) for h in y_labels]
            y_vals_err = [hists[h].GetBinError(ib, jb) for h in y_labels]
            sumw2 = [hists[h].GetSumw2()[hists[h].GetBin(ib, jb)] for h in y_labels]
            sumwx2 = [hists['%s_2' % h].GetBinContent(ib, jb) for h in y_labels]
            A = hists['A'].GetBinContent(ib, jb)
            B = hists['B'].GetBinContent(ib, jb)
            A_2 = hists['A_2'].GetBinContent(ib, jb)
            B_2 = hists['B_2'].GetBinContent(ib, jb)
            bin_str = '%.4g #leq %s < %.4g' % (hists[h].GetXaxis().GetBinLowEdge(ib), label_x, hists[h].GetXaxis().GetBinUpEdge(ib))
            bin_str +=', %.4g #leq %s < %.4g' % (hists[h].GetYaxis().GetBinLowEdge(jb), label_y, hists[h].GetYaxis().GetBinUpEdge(jb))
            bin_scalings[ib - 1].append(ParametrizeBin(x_vals, y_vals, y_vals_err, sumw2, A, B, A_2, B_2, '%s_%s_%i_%i' % (args.label, pm_label, ib - 1, jb - 1), makePlots=True, wsp=wsp, binStr=bin_str))
            scale = bin_scalings[ib - 1][jb - 1]
            for label, val in compute:
                scale_factor = (scale[0] + val * scale[1] + val * val * scale[2])
                hists[label].SetBinContent(ib, jb, hists[label].GetBinContent(ib, jb) * scale_factor)

    for ib in xrange(0, hists['nominal'].GetNbinsX()):
        for jb in xrange(0, hists['nominal'].GetNbinsY()):
            print '%-10i %-10i %-20.2f %15.2f +/- %-15.2f  %15.2f +/- %-15.2f' % (ib, jb, bin_scalings[ib][jb][0], bin_scalings[ib][jb][1], bin_scalings[ib][jb][3], bin_scalings[ib][jb][2], bin_scalings[ib][jb][4])

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

    legend.AddEntry(h_nominal, 'C_{3W} = 0.0', 'L')
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
    plot.DrawTitle(pads[0], 'W^{%s}#gamma^{}' %
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


for path, name, obj in hists.ListObjects():
    print (path, name, obj)
    WriteToTFile(obj, fout, path, name)
wsp.Write()
fout.Close()
