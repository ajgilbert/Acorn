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
parser.add_argument('--var', default='l0_pt', choices=['l0_pt', 'p0_chiso'])
parser.add_argument('--channel', default='m', choices=['e', 'm'])
parser.add_argument('--mc', action='store_true')

args = parser.parse_args()


fin = ROOT.TFile(args.input)
hists = TDirToNode(fin)
post = args.postfix

var = args.var
# var = 'p0_chiso'


def RebinHist(hist, bins):
    x = hist.Rebin(len(bins) - 1, "", array('d', bins))
    x.Copy(hist)


def ComputeSFs(node, strategy=0, var='l0j0_dphi', node_m=None):
    # Subtract all MC from the data - careful with GJets - assuming this is never
    # part of Total_J

    node['data_sub'] = node['data_obs'] - (node['Total_R'] + node['Total_E'] + node['Total_J'] - node['W_J'])

    # Strategy 0 - solve simultaneous
    if var == 'l0j0_dphi':
        obs_central = node['data_sub'].Integral(3, 28)
        obs_edge = node['data_sub'].Integral(1, 2) + node['data_sub'].Integral(29, 30)
        w_central = node['W_J'].Integral(3, 28)
        w_edge = node['W_J'].Integral(1, 2) + node['W_J'].Integral(29, 30)
        gjets_central = node['GJets_J'].Integral(3, 28)
        gjets_edge = node['GJets_J'].Integral(1, 2) + node['GJets_J'].Integral(29, 30)
    if var == 'puppi_met':
        obs_central = node['data_sub'].Integral(1, 7)
        obs_edge = node['data_sub'].Integral(8, 20)
        w_central = node['W_J'].Integral(1, 7)
        w_edge = node['W_J'].Integral(8, 20)
        gjets_central = node['GJets_J'].Integral(1, 7)
        gjets_edge = node['GJets_J'].Integral(8, 20)
    if strategy == 0:
        w_sf = (obs_central * (gjets_edge / gjets_central) - obs_edge) / (w_central * (gjets_edge / gjets_central) - w_edge)
        gjets_sf = (obs_edge - (w_edge * w_sf)) / gjets_edge
    elif strategy == 1:
        w_sf = 1.0
        gjets_sf = (obs_edge - (w_edge * 1.0)) / gjets_edge
    if strategy == 2:
        h_mu_sub = node_m['data_obs'] - (node_m['Total_R'] + node_m['Total_E'] + node_m['Total_J'] - node_m['W_J'])
        h_mu_w = node_m['W_J']
        w_sf = h_mu_sub.Integral() / h_mu_w.Integral()
        gjets_sf = (obs_edge - (w_edge * w_sf)) / gjets_edge
    return w_sf, gjets_sf


fout = ROOT.TFile(args.output, 'RECREATE')

res = {}

bins = [30, 40, 50, 60, 80, 100, 200]

for eb in ['barrel_%s' % args.channel, 'endcap_%s' % args.channel]:
    print '\n'
    for ib, b in enumerate(bins[:-1]):
        region_iso_t = '%s_iso_t_pt_%i_%i' % (eb, bins[ib], bins[ib + 1])
        region_iso_l = '%s_iso_l_pt_%i_%i' % (eb, bins[ib], bins[ib + 1])
        region_w_ctl = '%s_w_ctl_pt_%i_%i' % (eb, bins[ib], bins[ib + 1])
        print '>> %s, %s, %s' % (region_iso_t, region_iso_l, region_w_ctl)

        node_iso_t = hists[args.channel][region_iso_t]['l0met_mt']
        node_iso_l = hists[args.channel][region_iso_l]['l0met_mt']
        node_w_ctl = hists[args.channel][region_w_ctl]['l0met_mt']

        """
        Solve for (W+jets/Z+jets) and GJets SFs using the dphi distribution
        """
        var_w_ctl = 'l0j0_dphi'
        node_w_ctl_dphi = hists[args.channel][region_w_ctl][var_w_ctl]
        node_w_ctl_muons = hists['m'][region_w_ctl.replace('_e_', '_m_')][var_w_ctl]

        # Subtract all MC from the data - careful with GJets - assuming this is never
        # part of Total_J
        # node_w_ctl_dphi['data_sub'] = node_w_ctl_dphi['data_obs'] - (node_w_ctl_dphi['Total_R'] + node_w_ctl_dphi['Total_E'] + node_w_ctl_dphi['Total_J'] + node_w_ctl_dphi['GJets_J'])
        # w_ctl_dphi_central = node_w_ctl_dphi['data_sub'].Integral(3, 28)
        # w_ctl_dphi_edge = node_w_ctl_dphi['data_sub'].Integral(1, 2) + node_w_ctl_dphi['data_sub'].Integral(29, 30)
        # w_ctl_dphi_w_central = node_w_ctl_dphi['W_J'].Integral(3, 28)
        # w_ctl_dphi_w_edge = node_w_ctl_dphi['W_J'].Integral(1, 2) + node_w_ctl_dphi['W_J'].Integral(29, 30)
        # w_ctl_dphi_gjets_central = node_w_ctl_dphi['GJets_J'].Integral(3, 28)
        # w_ctl_dphi_gjets_edge = node_w_ctl_dphi['GJets_J'].Integral(1, 2) + node_w_ctl_dphi['GJets_J'].Integral(29, 30)

        # w_ctl_dphi_w_sf = ((w_ctl_dphi_w_central + w_ctl_dphi_central) / w_ctl_dphi_w_central)
        # w_ctl_dphi_edge = w_ctl_dphi_edge - (w_ctl_dphi_w_edge * (w_ctl_dphi_w_sf - 1.0))
        # w_ctl_dphi_gjets_sf = ((w_ctl_dphi_gjets_edge + w_ctl_dphi_edge) / w_ctl_dphi_gjets_edge)
        w_sf, gjets_sf = ComputeSFs(node_w_ctl_dphi, strategy=2, var=var_w_ctl, node_m=node_w_ctl_muons)
        print '>> dphi SFs: %f, %f' % (w_sf, gjets_sf)

        w_norm_w_ctl = node_w_ctl['W_J'].Integral()
        gjet_norm_w_ctl = node_w_ctl['GJets_J'].Integral()

        node_w_ctl['data_sub'] = node_w_ctl['data_obs'] - (node_w_ctl['Total_R'] + node_w_ctl['Total_E'] + node_w_ctl['Total_J'] + node_w_ctl['GJets_J'])
        data_sub_norm_w_ctl = node_w_ctl['data_sub'].Integral()
        # w_sf = ((w_norm_w_ctl + data_sub_norm_w_ctl) / w_norm_w_ctl)
        # gjet_sf = ((gjet_norm_w_ctl + data_sub_norm_w_ctl) / gjet_norm_w_ctl)
        # print w_sf, gjet_sf

        node_iso_t['Total_J'].Add(node_iso_t['W_J'], -1.0)
        node_iso_t['W_J'].Scale(w_sf)
        node_iso_t['GJets_J'].Scale(gjets_sf)
        node_iso_t['data_sub'] = node_iso_t['data_obs'] - (node_iso_t['Total_R'] + node_iso_t['Total_E'] + node_iso_t['Total_J'] + node_iso_t['W_J'] + node_iso_t['GJets_J'])
        node_iso_l['data_sub'] = node_iso_l['data_obs'] - (node_iso_l['Total_R'] + node_iso_l['Total_E'] + node_iso_l['Total_J'])
        data_sub_norm_iso_t = node_iso_t['data_sub'].Integral()
        data_sub_norm_iso_l = node_iso_l['data_sub'].Integral()

        print '>> Fake factor: %f' % (data_sub_norm_iso_t / data_sub_norm_iso_l)

        node_iso_t['data_sub'].Divide(node_iso_l['data_sub'])
        # node_iso_t['data_sub'].Print("range")
        # # Rebin if we're doing pT
        # if var == 'l0_pt':
        #     if 'barrel' in eb:
        #         # newbins = [30, 40, 60, 100, 200]
        #         newbins = [30, 40, 50, 60, 80, 100, 200]
        #     else:
        #         # newbins = [30, 40, 60, 100, 200]
        #         newbins = [30, 40, 50, 60, 80, 100, 200]
        #     node.ForEach(lambda x: RebinHist(x, newbins))

        # Subtract off the background
        # node['data_obs'].Print('range')
        # node['Total_R'].Print('range')
        # node['Total_E'].Print('range')
        # node['Total_J'].Print('range')
        # node['data_sub'] = node['data_obs'] - (node['Total_R'] + node['Total_E'] + node['Total_J'])

#     # Now take the ratio between regions
#     h_fr = hists[args.channel]['%s_iso_t%s' % (eb, post)][var]['data_sub'].Clone()
#     h_fr.Divide(hists[args.channel]['%s_iso_l%s' % (eb, post)][var]['data_sub'])
#     hint = None
#     if var == 'p0_chiso':
#         h_fr.Fit('pol1')
#         # func = ROOT.TF1('func', '[0]*TMath::Exp([1]*([2]-x))')
#         # h_fr.Fit('func')
#         hint = ROOT.TH1D("hint", "Fitted gaussian with .95 conf.band", 100, 0, 20)
#         ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hint, 0.68)
#         hint.SetStats(False)
#         hint.SetMarkerSize(0)
#         hint.SetFillColorAlpha(2, 0.2)
#         hint.Print('range')
#         # h_fr.GetFunction('pol2').Draw('E3SAME')

#     # Print the results
#     for ib in xrange(1, h_fr.GetNbinsX() + 1):
#         frac_err = 0.
#         if h_fr.GetBinContent(ib) > 0.:
#             frac_err = (h_fr.GetBinError(ib) / h_fr.GetBinContent(ib))
#         print '[%f,%f] %.3f %.3f %.3f' % (
#             h_fr.GetXaxis().GetBinLowEdge(ib), h_fr.GetXaxis().GetBinUpEdge(ib), h_fr.GetBinContent(ib), h_fr.GetBinError(ib), frac_err
#             )
#     if args.mc:
#         h_fr_mc = hists[args.channel]['%s_iso_l_sig_t' % eb][var]['W_J'].Clone()
#         h_fr_mc.Divide(hists[args.channel]['%s_iso_l_sig_l' % eb][var]['W_J'])
#         print '>> MC'
#         for ib in xrange(1, h_fr_mc.GetNbinsX() + 1):
#             frac_err = 0.
#             if h_fr_mc.GetBinContent(ib) > 0.:
#                 frac_err = (h_fr_mc.GetBinError(ib) / h_fr_mc.GetBinContent(ib))
#             print '[%f,%f] %.3f %.3f %.3f' % (
#                 h_fr_mc.GetXaxis().GetBinLowEdge(ib), h_fr_mc.GetXaxis().GetBinUpEdge(ib), h_fr_mc.GetBinContent(ib), h_fr_mc.GetBinError(ib), frac_err
#                 )
#     plot.Set(h_fr, MarkerSize=0.5)
#     canv = ROOT.TCanvas('lepton_fakes_%s_%s%s' % (args.year, eb, post), 'lepton_fakes_%s_%s%s' % (args.year, eb, post))
#     pads = plot.OnePad()
#     h_fr.Draw('E')
#     if hint is not None:
#         hint.Draw("E3SAME")
#         h_fr.Draw('ESAME')
#     if args.mc:
#         plot.Set(h_fr_mc, LineColor=2, MarkerColor=2, MarkerSize=0.5)
#         h_fr_mc.Draw('SAMEE')
#     h_fr.SetMaximum(10)
#     h_fr.SetMinimum(0)
#     h_fr.GetXaxis().SetTitle('Lepton p_{T} (GeV)')
#     if var == 'p0_chiso':
#         h_fr.GetXaxis().SetTitle('Photon I_{ch} (GeV)')
#     h_fr.GetYaxis().SetTitle('Fake ratio')
#     if args.mc:
#         legend = ROOT.TLegend(0.6, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')
#         legend.AddEntry(h_fr, 'Data', 'L')
#         legend.AddEntry(h_fr_mc, 'W+jets simulation', 'L')
#         legend.Draw()
#     canv.Print('.pdf')
#     canv.Print('.png')
#     fout.cd()
#     res[eb] = h_fr
#     ROOT.gDirectory.WriteTObject(h_fr, eb)


# h2d = ROOT.TH2F('lepton_fakes', '', 34, 30, 200, 2, array('d', [0., 1.4442, 2.5]))

# for i in xrange(1, h2d.GetNbinsX() + 1):
#     x = h2d.GetXaxis().GetBinCenter(i)
#     h2d.SetBinContent(i, 1, res['barrel_%s' % args.channel].GetBinContent(res['barrel_%s' % args.channel].GetXaxis().FindFixBin(x)))
#     h2d.SetBinContent(i, 2, res['endcap_%s' % args.channel].GetBinContent(res['endcap_%s' % args.channel].GetXaxis().FindFixBin(x)))
#     h2d.SetBinError(i, 1, res['barrel_%s' % args.channel].GetBinError(res['barrel_%s' % args.channel].GetXaxis().FindFixBin(x)))
#     h2d.SetBinError(i, 2, res['endcap_%s' % args.channel].GetBinError(res['endcap_%s' % args.channel].GetXaxis().FindFixBin(x)))

# h2d_syst_extrap = h2d.Clone()
# h2d_stat_syst = h2d.Clone()

# ROOT.gDirectory.WriteTObject(h2d, 'lepton_fakes')


fout.Close()