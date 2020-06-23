import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

import CombineHarvester.CombineTools.plotting as plot
# import json
# import sys
import math
from pprint import pprint
from collections import defaultdict
import argparse
from array import array
from Acorn.Analysis.analysis import *

ROOT.TH1.SetDefaultSumw2(True)
plot.ModTDRStyle(height=300)

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', default='')
parser.add_argument('-o', '--output', default='output_2016_photon_fakes_ratios_m.root')
parser.add_argument('--year', default='2016')
parser.add_argument('--postfix', default='')
parser.add_argument('--var', default='p0_pt', choices=['p0_pt', 'p0_chiso'])
parser.add_argument('--channel', default='m', choices=['e', 'm'])
parser.add_argument('--mc', action='store_true')

args = parser.parse_args()


fin = ROOT.TFile(args.input)
hists = TDirToNode(fin)
post = args.postfix

var = args.var
# var = 'p0_chiso'
add_channel = 'e'

def RebinHist(hist, bins):
    x = hist.Rebin(len(bins) - 1, "", array('d', bins))
    x.Copy(hist)


fout = ROOT.TFile(args.output, 'RECREATE')

res = {}


def DoMCRatio(num, den, rebin=None):
    num.Print('range')
    den.Print('range')
    h_fr_mc = num.Clone()
    h_den = den.Clone()
    if rebin is not None:
        h_fr_mc.Rebin(rebin)
        h_den.Rebin(rebin)
    h_fr_mc.Divide(h_den)
    # print '>> MC'
    for ib in xrange(1, h_fr_mc.GetNbinsX() + 1):
        frac_err = 0.
        if h_fr_mc.GetBinContent(ib) > 0.:
            frac_err = (h_fr_mc.GetBinError(ib) / h_fr_mc.GetBinContent(ib))
            print '[%f,%f] %.3f %.3f %.3f' % (
                h_fr_mc.GetXaxis().GetBinLowEdge(ib), h_fr_mc.GetXaxis().GetBinUpEdge(ib), h_fr_mc.GetBinContent(ib), h_fr_mc.GetBinError(ib), frac_err
                )
    return h_fr_mc


# for eb in ['barrel_%s' % args.channel, 'barrel1_%s' % args.channel, 'barrel2_%s' % args.channel, 'endcap_%s' % args.channel, 'endcap1_%s' % args.channel, 'endcap2_%s' % args.channel]:
for eb in ['barrel_%s' % args.channel, 'endcap_%s' % args.channel]:
    for sel in ['%s_iso_t_sig_t%s' % (eb, post), '%s_iso_t_sig_l%s' % (eb, post), '%s_iso_l_sig_l%s' % (eb, post), '%s_iso_l_sig_t%s' % (eb, post), '%s_sig_l%s' % (eb, post), '%s_sig_t%s' % (eb, post)]:
        node = hists[args.channel][sel][var]
        sel_alt = sel.replace(eb, eb.replace('_%s' % args.channel, '_%s' % add_channel))
        print sel_alt
        node_alt = hists[add_channel][sel_alt][var]

        # Rebin if we're doing pT
        if var == 'p0_pt':
            if 'barrel_' in eb or 'barrel1_' in eb:
                # newbins = [30, 40, 60, 100, 200]
                newbins = [30, 40, 50, 60, 100, 150, 300]
            elif 'barrel2_' in eb:
                newbins = [30, 40, 50, 60, 100, 300]
            elif 'endcap_' in eb or 'endcap1_' in eb:
                # newbins = [30, 40, 60, 100, 200]
                newbins = [30, 40, 60, 100, 300]
            elif 'endcap2_' in eb:
                newbins = [30, 40, 60, 300]
            node.ForEach(lambda x: RebinHist(x, newbins))
            if add_channel is not None:
                node_alt.ForEach(lambda x: RebinHist(x, newbins))

        # Subtract off the background
        node['data_sub'] = node['data_obs'] - (node['Total_R'] + node['Total_E'])
        if add_channel is not None: 
            node['data_sub'] += (node_alt['data_obs'] - (node_alt['Total_R'] + node_alt['Total_E']))

    # Now take the ratio between regions
    h_fr = hists[args.channel]['%s_sig_t%s' % (eb, post)][var]['data_sub'].Clone()
    h_fr.Divide(hists[args.channel]['%s_sig_l%s' % (eb, post)][var]['data_sub'])
    hint = None
    if var == 'p0_chiso':
        h_fr.Fit('pol1')
        # func = ROOT.TF1('func', '[0]*TMath::Exp([1]*([2]-x))')
        # h_fr.Fit('func')
        hint = ROOT.TH1D("hint", "Fitted gaussian with .95 conf.band", 100, 0, 20)
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hint, 0.68)
        hint.SetStats(False)
        hint.SetMarkerSize(0)
        hint.SetFillColorAlpha(2, 0.2)
        # hint.Print('range')
        # h_fr.GetFunction('pol2').Draw('E3SAME')

    # Print the results
    for ib in xrange(1, h_fr.GetNbinsX() + 1):
        frac_err = 0.
        if h_fr.GetBinContent(ib) > 0.:
            frac_err = (h_fr.GetBinError(ib) / h_fr.GetBinContent(ib))
        print '[%f,%f] %.3f %.3f %.3f' % (
            h_fr.GetXaxis().GetBinLowEdge(ib), h_fr.GetXaxis().GetBinUpEdge(ib), h_fr.GetBinContent(ib), h_fr.GetBinError(ib), frac_err
            )
    if args.mc:
        print '>> MC'
        h_fr_mc = DoMCRatio(hists[args.channel]['%s_iso_l_sig_t%s' % (eb, post)][var]['W_J'], hists[args.channel]['%s_iso_l_sig_l%s' % (eb, post)][var]['W_J'], rebin=1)
        h_fr_mc_all = DoMCRatio(hists[args.channel]['%s_sig_t%s' % (eb, post)][var]['W_J'], hists[args.channel]['%s_sig_l%s' % (eb, post)][var]['W_J'], rebin=2)
        print '>> MC (h)'
        h_fr_mc_h = DoMCRatio(hists[args.channel]['%s_iso_l_sig_t%s' % (eb, post)][var]['W_H'], hists[args.channel]['%s_iso_l_sig_l%s' % (eb, post)][var]['W_H'])
        print '>> MC (truth)'
        h_fr_mc_truth = DoMCRatio(hists[args.channel]['%s_iso_t_sig_t%s' % (eb, post)][var]['W_J'], hists[args.channel]['%s_iso_t_sig_l%s' % (eb, post)][var]['W_J'])
        h_fr_mc_truth_all = DoMCRatio(hists[args.channel]['%s_sig_t%s' % (eb, post)][var]['Total_J'], hists[args.channel]['%s_sig_l%s' % (eb, post)][var]['Total_J'], rebin=1)
        print '>> MC (h - truth)'
        h_fr_mc_h_truth = DoMCRatio(hists[args.channel]['%s_sig_t%s' % (eb, post)][var]['W_H'], hists[args.channel]['%s_sig_l%s' % (eb, post)][var]['W_H'], rebin=2)
        print '>> MC (p - truth)'
        h_fr_mc_p_truth = DoMCRatio(hists[args.channel]['%s_sig_t%s' % (eb, post)][var]['W_P'], hists[args.channel]['%s_sig_l%s' % (eb, post)][var]['W_P'], rebin=2)
        print '>> MC (m - truth)'
        h_fr_mc_m_truth = DoMCRatio(hists[args.channel]['%s_sig_t%s' % (eb, post)][var]['W_M'], hists[args.channel]['%s_sig_l%s' % (eb, post)][var]['W_M'], rebin=2)
    plot.Set(h_fr, MarkerSize=0.5)
    canv = ROOT.TCanvas('photon_fakes_%s_%s%s' % (args.year, eb, post), 'photon_fakes_%s_%s%s' % (args.year, eb, post))
    pads = plot.OnePad()
    h_fr.Draw('E')
    if hint is not None:
        hint.Draw("E3SAME")
        h_fr.Draw('ESAME')
    if args.mc:
        plot.Set(h_fr_mc, LineColor=2, MarkerColor=2, MarkerSize=0.5)
        plot.Set(h_fr_mc_all, LineColor=2, MarkerColor=2, MarkerSize=0.5)
        plot.Set(h_fr_mc_truth, LineColor=4, MarkerColor=4, MarkerSize=0.5)
        plot.Set(h_fr_mc_truth_all, LineColor=4, MarkerColor=4, MarkerSize=0.5)
        plot.Set(h_fr_mc_h, LineColor=30, MarkerColor=30, MarkerSize=0.5)
        plot.Set(h_fr_mc_h_truth, LineColor=6, MarkerColor=6, MarkerSize=0.5)
        plot.Set(h_fr_mc_p_truth, LineColor=8, MarkerColor=8, MarkerSize=0.5)
        plot.Set(h_fr_mc_m_truth, LineColor=28, MarkerColor=28, MarkerSize=0.5)

        # hint_mc = None
        # if var == 'p0_chiso':
        #     h_fr_mc_truth_all.Fit('pol1')
        #     # func = ROOT.TF1('func', '[0]*TMath::Exp([1]*([2]-x))')
        #     # h_fr.Fit('func')
        #     hint_mc = ROOT.TH1D("hint_mc", "Fitted gaussian with .95 conf.band", 100, 0, 20)
        #     ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hint_mc, 0.68)
        #     hint_mc.SetStats(False)
        #     hint_mc.SetMarkerSize(0)
        #     hint_mc.SetFillColorAlpha(4, 0.2)
        #     hint_mc.SetLineColor(4)
        #     # hint_mc.Print('range')
        # if hint_mc is not None:
        #     hint_mc.Draw("E3SAME")
        #     h_fr.Draw('ESAME')

        h_fr_mc.Draw('SAMEE')
        # h_fr_mc_truth.Draw('SAMEE')
        # h_fr_mc_all.Draw('SAMEE')
        # h_fr_mc_truth_all.Rebin(2)
        # h_fr_mc_truth_all.Scale(0.5)
        h_fr_mc_truth_all.Draw('SAMEE')
        # h_fr_mc_h.Draw('SAMEE')
        # h_fr_mc_h_truth.Draw('SAMEE')
        # h_fr_mc_p_truth.Draw('SAMEE')
        # h_fr_mc_m_truth.Draw('SAMEE')
    h_fr.SetMaximum(2.1)
    h_fr.SetMinimum(0.1)
    h_fr.GetXaxis().SetTitle('Photon p_{T} (GeV)')
    if var == 'p0_chiso':
        h_fr.GetXaxis().SetTitle('Photon I_{ch} (GeV)')
    h_fr.GetYaxis().SetTitle('Fake ratio')
    if args.mc:
        legend = ROOT.TLegend(0.6, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')
        legend.AddEntry(h_fr, 'Data', 'L')
        # legend.AddEntry(h_fr_mc, 'W+jets (anti-iso)', 'L')
        # legend.AddEntry(h_fr_mc_truth, 'W+jets (iso)', 'L')
        # legend.AddEntry(h_fr_mc_all, 'W+jets (anti-iso)', 'L')
        legend.AddEntry(h_fr_mc_truth_all, 'W+jets (iso)', 'L')
        # legend.AddEntry(h_fr_mc_h, 'W+jets [H] (anti-iso)', 'L')
        # legend.AddEntry(h_fr_mc_h_truth, 'W+jets [H] (iso)', 'L')
        # legend.AddEntry(h_fr_mc_p_truth, 'W+jets [P] (iso)', 'L')
        # legend.AddEntry(h_fr_mc_m_truth, 'W+jets [M] (iso)', 'L')
        legend.Draw()
    canv.Print('.pdf')
    canv.Print('.png')
    fout.cd()
    res[eb] = h_fr
    ROOT.gDirectory.WriteTObject(h_fr, eb)

use4bins = True

if use4bins:
    h2d = ROOT.TH2F('photon_fakes', '', 54, 30, 300, 4, array('d', [0., 1.0, 1.4442, 2.1, 2.5]))
else:
    h2d = ROOT.TH2F('photon_fakes', '', 54, 30, 300, 2, array('d', [0., 1.4442, 2.5]))

h2d_index = h2d.Clone()
curr_index = 0
for j in xrange(1, h2d.GetNbinsY() + 1):
    curr_index += 1
    if j == 1:
        curr_hist = res['barrel1_%s' % args.channel]
    if j == 2:
        curr_hist = res['barrel2_%s' % args.channel]
    if j == 3:
        curr_hist = res['endcap1_%s' % args.channel]
    if j == 4:
        curr_hist = res['endcap2_%s' % args.channel]
    curr_x_bin = 1
    for i in xrange(1, h2d.GetNbinsX() + 1):
        x = h2d.GetXaxis().GetBinCenter(i)
        xbin = curr_hist.GetXaxis().FindFixBin(x)
        if xbin > curr_x_bin:
            curr_x_bin = xbin
            curr_index += 1
        h2d_index.SetBinContent(i, j, float(curr_index))


for i in xrange(1, h2d.GetNbinsX() + 1):
    x = h2d.GetXaxis().GetBinCenter(i)
    if use4bins:
        h2d.SetBinContent(i, 1, res['barrel1_%s' % args.channel].GetBinContent(res['barrel1_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinContent(i, 2, res['barrel2_%s' % args.channel].GetBinContent(res['barrel2_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinContent(i, 3, res['endcap1_%s' % args.channel].GetBinContent(res['endcap1_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinContent(i, 4, res['endcap2_%s' % args.channel].GetBinContent(res['endcap2_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinError(i, 1, res['barrel1_%s' % args.channel].GetBinError(res['barrel1_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinError(i, 2, res['barrel2_%s' % args.channel].GetBinError(res['barrel2_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinError(i, 3, res['endcap1_%s' % args.channel].GetBinError(res['endcap1_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinError(i, 4, res['endcap2_%s' % args.channel].GetBinError(res['endcap2_%s' % args.channel].GetXaxis().FindFixBin(x)))
    else:
        h2d.SetBinContent(i, 1, res['barrel_%s' % args.channel].GetBinContent(res['barrel_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinContent(i, 2, res['endcap_%s' % args.channel].GetBinContent(res['endcap_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinError(i, 1, res['barrel_%s' % args.channel].GetBinError(res['barrel_%s' % args.channel].GetXaxis().FindFixBin(x)))
        h2d.SetBinError(i, 2, res['endcap_%s' % args.channel].GetBinError(res['endcap_%s' % args.channel].GetXaxis().FindFixBin(x)))


"""
These are the linear extrapolations in charged iso to 0.5 GeV:
BARREL - 2016
30-40:
 fSumw[3]=1.2113, x=0.5, error=0.0705191
40-60:
 fSumw[3]=0.995309, x=0.5, error=0.0586166
60-100:
 fSumw[3]=0.697634, x=0.5, error=0.0840808
100-200:
 fSumw[3]=0.280704, x=0.5, error=0.0550755

[30.000000,40.000000] 1.004965 0.018217 0.018127
[40.000000,60.000000] 0.887880 0.017478 0.019685
[60.000000,100.000000] 0.657777 0.015734 0.023920
[100.000000,200.000000] 0.363980 0.015708 0.043155
Rel uncerts (+ stat err):
 - 1.205315608
 - 1.120994954
 - 1.060593484
 - 0.771207209

ENDCAP - 2016
30-40:
 fSumw[3]=0.771773, x=0.5, error=0.0696247
40-60:
 fSumw[3]=0.703544, x=0.5, error=0.0832643
60-100:
 fSumw[3]=0.898874, x=0.5, error=0.100833
 100-200:
 fSumw[3]=0.703094, x=0.5, error=0.222491

[30.000000,40.000000] 0.613333 0.015178 0.024747
[40.000000,60.000000] 0.603833 0.016708 0.027670
[60.000000,100.000000] 0.679418 0.024028 0.035366
[100.000000,200.000000] 0.833903 0.056843 0.068165

Rel uncerts (+ stat err):
 - 1.258326227
 - 1.165130094
 - 1.323005867
 - 0.843136432


BARREL - 2017
30-40:
 fSumw[3]=1.15364, x=0.5, error=0.055392
40-60:
 fSumw[3]=0.875979, x=0.5, error=0.0489659
60-100:
 fSumw[3]=0.677147, x=0.5, error=0.0544896
100-200:
 fSumw[3]=0.45907, x=0.5, error=0.0802994

[30.000000,40.000000] 0.942330 0.013639 0.014473
[40.000000,60.000000] 0.827754 0.013030 0.015742
[60.000000,100.000000] 0.644211 0.012020 0.018658
[100.000000,200.000000] 0.477286 0.015821 0.033149
Rel uncerts (+ stat err):
 - 1.220379273
 - 1.058260063
 - 1.051126106
 - 0.961834204

ENDCAP - 2017
30-40:
 fSumw[3]=0.805178, x=0.5, error=0.0485607
40-60:
 fSumw[3]=0.886195, x=0.5, error=0.0707819
60-100:
 fSumw[3]=1.01639, x=0.5, error=0.115265
 100-200:
 fSumw[5]=1.07788, x=0.9, error=0.301499

[30.000000,40.000000] 0.591468 0.014958 0.025290
[40.000000,60.000000] 0.690967 0.019344 0.027996
[60.000000,100.000000] 0.864095 0.030710 0.035540
[100.000000,200.000000] 1.260906 0.088633 0.070293

Rel uncerts (+ stat err):
 - 1.361321323
 - 1.282543161
 - 1.176247982
 - 0.854845643


BARREL - 2018
30-40:
 fSumw[3]=1.21661, x=0.5, error=0.0256668
40-60:
 fSumw[3]=0.96831, x=0.5, error=0.0325356
60-100:
 fSumw[3]=0.63544, x=0.5, error=0.0410728
100-200:
 fSumw[3]=0.502062, x=0.5, error=0.0572004

[30.000000,40.000000] 1.012072 0.014289 0.014119
[40.000000,60.000000] 0.875597 0.012068 0.013782
[60.000000,100.000000] 0.663844 0.010914 0.016440
[100.000000,200.000000] 0.474364 0.013239 0.027910
Rel uncerts (+ stat err):
 - 1.20209827
 - 1.10588547
 - 0.957212839
 - 1.05838976

ENDCAP - 2018
30-40:
 fSumw[3]=0.828486, x=0.5, error=0.0500235
40-60:
 fSumw[3]=0.846882, x=0.5, error=0.0597682
60-100:
 fSumw[3]=0.807961, x=0.5, error=0.0829876
100-200:
 fSumw[3]=1.01117, x=0.5, error=0.201686

[30.000000,40.000000] 0.555193 0.011527 0.020761
[40.000000,60.000000] 0.585242 0.013790 0.023563
[60.000000,100.000000] 0.788115 0.023815 0.030217
[100.000000,200.000000] 1.201661 0.070848 0.058959

Rel uncerts (+ stat err):
 - 1.492248641
 - 1.447062924
 - 1.025181604
 - 0.841476922


Overall:
Barrel:
 - 30-40:  1.20
 - 40-60:  1.10
 - 60-100: 1.05
 - 100+:   1.20
 Endcap:
  - 30-40:  1.50
  - 40-60:  1.45
  - 60-100: 1.30
  - 100+:   1.15
"""
syst_vals_barrel = [0.20, 0.10, 0.05, 0.20]
syst_vals_endcap = [0.50, 0.45, 0.30, 0.15]
h_syst_barrel = ROOT.TH1F('h_syst_barrel', '', 4, array('d', [30, 40, 60, 100, 300]))
h_syst_endcap = ROOT.TH1F('h_syst_endcap', '', 4, array('d', [30, 40, 60, 100, 300]))

for ib in xrange(1, h_syst_barrel.GetNbinsX() + 1):
    h_syst_barrel.SetBinContent(ib, syst_vals_barrel[ib - 1])
    h_syst_endcap.SetBinContent(ib, syst_vals_endcap[ib - 1])

h2d_syst_extrap = h2d.Clone()
h2d_stat_syst = h2d.Clone()

for i in xrange(1, h2d_syst_extrap.GetNbinsX() + 1):
    x = h2d_syst_extrap.GetXaxis().GetBinCenter(i)
    if use4bins:
        h2d_syst_extrap.SetBinError(i, 1, h2d_syst_extrap.GetBinContent(i, 1) * h_syst_barrel.GetBinContent(h_syst_barrel.GetXaxis().FindFixBin(x)))
        h2d_syst_extrap.SetBinError(i, 2, h2d_syst_extrap.GetBinContent(i, 2) * h_syst_barrel.GetBinContent(h_syst_barrel.GetXaxis().FindFixBin(x)))
        h2d_syst_extrap.SetBinError(i, 3, h2d_syst_extrap.GetBinContent(i, 3) * h_syst_endcap.GetBinContent(h_syst_endcap.GetXaxis().FindFixBin(x)))
        h2d_syst_extrap.SetBinError(i, 4, h2d_syst_extrap.GetBinContent(i, 4) * h_syst_endcap.GetBinContent(h_syst_endcap.GetXaxis().FindFixBin(x)))
        h2d_stat_syst.SetBinError(i, 1, math.sqrt(math.pow(h2d.GetBinError(i, 1), 2) + math.pow(h2d_syst_extrap.GetBinError(i, 1), 2)))
        h2d_stat_syst.SetBinError(i, 2, math.sqrt(math.pow(h2d.GetBinError(i, 2), 2) + math.pow(h2d_syst_extrap.GetBinError(i, 2), 2)))
        h2d_stat_syst.SetBinError(i, 3, math.sqrt(math.pow(h2d.GetBinError(i, 3), 2) + math.pow(h2d_syst_extrap.GetBinError(i, 3), 2)))
        h2d_stat_syst.SetBinError(i, 4, math.sqrt(math.pow(h2d.GetBinError(i, 4), 2) + math.pow(h2d_syst_extrap.GetBinError(i, 4), 2)))
    else:
        h2d_syst_extrap.SetBinError(i, 1, h2d_syst_extrap.GetBinContent(i, 1) * h_syst_barrel.GetBinContent(h_syst_barrel.GetXaxis().FindFixBin(x)))
        h2d_syst_extrap.SetBinError(i, 2, h2d_syst_extrap.GetBinContent(i, 2) * h_syst_endcap.GetBinContent(h_syst_endcap.GetXaxis().FindFixBin(x)))
        h2d_stat_syst.SetBinError(i, 1, math.sqrt(math.pow(h2d.GetBinError(i, 1), 2) + math.pow(h2d_syst_extrap.GetBinError(i, 1), 2)))
        h2d_stat_syst.SetBinError(i, 2, math.sqrt(math.pow(h2d.GetBinError(i, 2), 2) + math.pow(h2d_syst_extrap.GetBinError(i, 2), 2)))

h2d.Print('range')
h2d_syst_extrap.Print('range')
h2d_stat_syst.Print('range')

ROOT.gDirectory.WriteTObject(h2d, 'photon_fakes')
ROOT.gDirectory.WriteTObject(h2d_index, 'photon_fakes_index')
ROOT.gDirectory.WriteTObject(h2d_syst_extrap, 'photon_fakes_syst')
ROOT.gDirectory.WriteTObject(h2d_stat_syst, 'photon_fakes_stat_syst')

fout.Close()