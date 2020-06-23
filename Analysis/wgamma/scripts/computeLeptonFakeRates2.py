import argparse
from array import array
from Acorn.Analysis.analysis import *
import CombineHarvester.CombineTools.plotting as plot
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# import math
# from pprint import pprint
# from collections import defaultdict

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
parser.add_argument('--syst', type=int, default=0)

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
        l_max = 3
        h_max = 27
        obs_central = node['data_sub'].Integral(l_max + 1, h_max - 1)
        obs_edge = node['data_sub'].Integral(1, l_max) + node['data_sub'].Integral(h_max, 30)
        w_central = node['W_J'].Integral(l_max + 1, h_max - 1)
        w_edge = node['W_J'].Integral(1, l_max) + node['W_J'].Integral(h_max, 30)
        gjets_central = node['GJets_J'].Integral(l_max + 1, h_max - 1)
        gjets_edge = node['GJets_J'].Integral(1, l_max) + node['GJets_J'].Integral(h_max, 30)
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

bins = [30, 50, 100]

x_bins = [30, 40, 50, 60, 70, 80, 90, 100]
z_bins = [0, 10, 20, 30, 40, 50, 60, 70, 80]
if args.channel == 'e':
    y_bins = [0., 1.4442, 2.5]
else:
    y_bins = [0., 1.5, 2.5]
h2d = ROOT.TH3F('lepton_fakes', '', len(x_bins) - 1, array('d', x_bins), len(y_bins) - 1, array('d', y_bins), len(z_bins) - 1, array('d', z_bins))

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
        # node_w_ctl_muons = hists[args.channel][region_w_ctl][var_w_ctl]

        w_sf, gjets_sf = ComputeSFs(node_w_ctl_dphi, strategy=2, var=var_w_ctl, node_m=node_w_ctl_muons)
        if args.channel == 'm':
            gjets_sf = 1.0
        gjets_sf = 1.0

        if args.syst == 1:
            w_sf = ((w_sf - 1.0) / 2.) + 1.0
            # gjets_sf = ((gjets_sf - 1.0) / 2.) + 1.0
        if args.syst == 2:
            w_sf = ((w_sf - 1.0) / 2.) + w_sf
            # gjets_sf = ((gjets_sf - 1.0) / 2.) + gjets_sf

        print '>> dphi SFs: %.2f, %.2f' % (w_sf, gjets_sf)

        w_norm_w_ctl = node_w_ctl['W_J'].Integral()
        gjet_norm_w_ctl = node_w_ctl['GJets_J'].Integral()

        node_w_ctl['data_sub'] = node_w_ctl['data_obs'] - (node_w_ctl['Total_R'] + node_w_ctl['Total_E'] + node_w_ctl['Total_J'] + node_w_ctl['GJets_J'])
        data_sub_norm_w_ctl = node_w_ctl['data_sub'].Integral()
        node_iso_t['Total_J'].Add(node_iso_t['W_J'], -1.0)
        node_iso_t['W_J'].Scale(w_sf)
        node_iso_t['GJets_J'].Scale(gjets_sf)
        node_iso_t['data_sub'] = node_iso_t['data_obs'] - (node_iso_t['Total_R'] + node_iso_t['Total_E'] + node_iso_t['Total_J'] + node_iso_t['W_J'] + node_iso_t['GJets_J'])
        node_iso_l['data_sub'] = node_iso_l['data_obs'] - (node_iso_l['Total_R'] + node_iso_l['Total_E'] + node_iso_l['Total_J'])

        h_ratio = node_iso_t['data_sub'].Clone()
        h_ratio.Divide(node_iso_l['data_sub'])
        h_ratio.Print("range")

        h_data_sub_norm_iso_t = VariableRebin(node_iso_t['data_sub'], [node_iso_t['data_sub'].GetXaxis().GetXmin(), node_iso_t['data_sub'].GetXaxis().GetXmax()])
        h_data_sub_norm_iso_l = VariableRebin(node_iso_l['data_sub'], [node_iso_l['data_sub'].GetXaxis().GetXmin(), node_iso_l['data_sub'].GetXaxis().GetXmax()])
        h_data_sub_norm_iso_t.Divide(h_data_sub_norm_iso_l)
        fake_factor = h_data_sub_norm_iso_t.GetBinContent(1)
        fake_factor_err = h_data_sub_norm_iso_t.GetBinError(1)

        print '>> Fake factor: %.2f +/- %.2f' % (fake_factor, fake_factor_err)

        iy = 1 if 'barrel' in eb else 2
        for ix in xrange(1, h2d.GetNbinsX() + 1):
            for iz in xrange(1, h2d.GetNbinsZ() + 1):
                xcenter = h2d.GetXaxis().GetBinCenter(ix)
                if xcenter > bins[ib] and xcenter < bins[ib + 1]:
                    if args.channel == 'e':
                        # take the FFs from h_ratio vs mT
                        h_ratio_bin = h_ratio.GetXaxis().FindFixBin(h2d.GetZaxis().GetBinCenter(iz))
                        h2d.SetBinContent(ix, iy, iz, h_ratio.GetBinContent(h_ratio_bin))
                    else:
                        h2d.SetBinContent(ix, iy, iz, fake_factor)
        # node_iso_t['data_sub'].Divide(node_iso_l['data_sub'])


h2d.Print('range')

plot_name = 'lepton_fakes' if args.syst == 0 else 'lepton_fakes_%i' % args.syst

ROOT.gDirectory.WriteTObject(h2d, plot_name)


fout.Close()
