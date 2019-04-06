#!/usr/bin/env python
import ROOT
import imp
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.workspaceTools import *
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()

fin = ROOT.TFile('hvm_corrections_2018_v2.root')
w = fin.Get('w')

def Compare(w, output, fn1, fn2, bins, var, other_vars= {}, line_pos=None):
    label_vars = []
    for key, val in other_vars.iteritems():
        w.var(key).setVal(val)
        label_vars.append('%s=%g' % (key, val))
    h1 = w.function(fn1).createHistogram(fn1, w.var(var),
            ROOT.RooFit.Binning(*bins),
            ROOT.RooFit.Scaling(False)
        )
    h2 = w.function(fn2).createHistogram(fn2, w.var(var),
            ROOT.RooFit.Binning(*bins),
            ROOT.RooFit.Scaling(False)
        )
    canv = ROOT.TCanvas(output, output)
    pads = plot.TwoPadSplit(0.30, 0.01, 0.01)
    pads[0].cd()
    pads[0].SetGrid(1, 1)
    plot.Set(h1, LineColor=ROOT.kBlack, LineWidth=2)
    plot.Set(h1.GetYaxis(), Title='Efficiency')
    plot.Set(h2, LineColor=ROOT.kRed, LineWidth=2)

    for i in xrange(1, h1.GetNbinsX()+1):
        h1.SetBinError(i, 0.)
    for i in xrange(1, h2.GetNbinsX()+1):
        h2.SetBinError(i, 0.)
    h1.Draw('L')
    h2.Draw('LSAME')
    ratio = h1.Clone()
    ratio.Divide(h2)


    legend = ROOT.TLegend(0.18, 0.82, 0.6, 0.93, '', 'NBNDC')
    legend.AddEntry(h1, fn1, 'L')
    legend.AddEntry(h2, fn2, 'L')
    legend.Draw()
    print plot.GetPadYMax(pads[0])
    plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.25)
    plot.DrawTitle(pads[0], ','.join(label_vars), 1)

    line = ROOT.TLine()
    plot.Set(line, LineColor=12, LineStyle=4, LineWidth=2)
    if line_pos is not None:
        plot.DrawVerticalLine(pads[0], line, line_pos)

    pads[1].cd()
    pads[1].SetGrid(1, 1)
    ratio.Draw('L')
    plot.SetupTwoPadSplitAsRatio(
        pads, plot.GetAxisHist(
            pads[0]), plot.GetAxisHist(pads[1]), 'Ratio', True, 0.91, 1.09)
    if line_pos is not None:
        plot.DrawVerticalLine(pads[1], line, line_pos)
    canv.Print('.pdf')
    canv.Print('.png')

w.var('rho_eta').setVal(1.3)

Compare(w, 'compare_data_mc_rho_iso_eff_inc', 'rho_iso_data_eff_etainc', 'rho_iso_mc_eff_etainc', bins=[160, 20, 100], var='rho_pt', other_vars={'rho_eta':1.3}, line_pos=20.)
Compare(w, 'compare_SF_inc_eta1', 'rhoiso_ratio_etainc', 'rhoiso_ratio', bins=[160, 20, 100], var='rho_pt', other_vars={'rho_eta':0.5}, line_pos=20.)
Compare(w, 'compare_SF_inc_eta1p56', 'rhoiso_ratio_etainc', 'rhoiso_ratio', bins=[160, 20, 100], var='rho_pt', other_vars={'rho_eta':1.2}, line_pos=20.)
Compare(w, 'compare_SF_inc_eta2p1', 'rhoiso_ratio_etainc', 'rhoiso_ratio', bins=[160, 20, 100], var='rho_pt', other_vars={'rho_eta':2.}, line_pos=20.)
Compare(w, 'compare_SF_inc_eta2p5', 'rhoiso_ratio_etainc', 'rhoiso_ratio', bins=[160, 20, 100], var='rho_pt', other_vars={'rho_eta':2.3}, line_pos=20.)
