import ROOT
import CombineHarvester.CombineTools.plotting as plot
import json
# import sys
from pprint import pprint
from collections import defaultdict
import argparse
from array import array
from Acorn.Analysis.analysis import *

ROOT.TH1.SetDefaultSumw2(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()

parser = argparse.ArgumentParser()

parser.add_argument('--year', default='2016', choices=['2016', '2017', '2018'])
parser.add_argument('--channel', default='m', choices=['e', 'm'])

args = parser.parse_args()


fin = ROOT.TFile('output_%s_photon_fakes.root' % args.year)
hists = TDirToNode(fin)

var = 'p0_pt'
# var = 'p0_chiso'

def RebinHist(hist, bins):
    x = hist.Rebin(len(bins) - 1, "", array('d', bins))
    x.Copy(hist)


fout = ROOT.TFile('output_%s_photon_fakes_ratios_%s.root' % (args.year, args.channel), 'RECREATE')

res = {}

for eb in ['barrel_%s' % args.channel, 'endcap_%s' % args.channel]:
    for sel in ['%s_iso_t_sig_l' % eb, '%s_iso_l_sig_l' % eb, '%s_iso_l_sig_t' % eb]:
        node = hists[args.channel][sel][var]
        if 'barrel' in eb:
            newbins = [30, 35, 40, 50, 60, 80, 100, 300]
        else:
            newbins = [30, 40, 60, 300]

        # newbins = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        node.ForEach(lambda x: RebinHist(x, newbins))
        node['data_sub'] = node['data_obs'] - (node['Total_R'] + node['Total_E'])
    h_fr = hists[args.channel]['%s_iso_l_sig_t' % eb][var]['data_sub'].Clone()
    h_fr.Divide(hists[args.channel]['%s_iso_l_sig_l' % eb][var]['data_sub'])

    h_fr_mc = hists[args.channel]['%s_iso_l_sig_t' % eb][var]['W_F'].Clone()
    h_fr_mc.Divide(hists[args.channel]['%s_iso_l_sig_l' % eb][var]['W_F'])

    plot.Set(h_fr, MarkerSize=0.5)
    canv = ROOT.TCanvas('photon_fakes_%s_%s' % (args.year, eb), 'photon_fakes_%s_%s' % (args.year, eb))
    pads = plot.OnePad()
    h_fr.Draw('E')
    plot.Set(h_fr_mc, LineColor=2, MarkerColor=2)
    h_fr_mc.Draw('SAMEE')
    h_fr.SetMaximum(1.2)
    h_fr.SetMinimum(0.1)
    h_fr.GetXaxis().SetTitle('Photon p_{T} (GeV)')
    h_fr.GetYaxis().SetTitle('Fake ratio')
    canv.Print('.pdf')
    canv.Print('.png')
    fout.cd()
    res[eb] = h_fr
    ROOT.gDirectory.WriteTObject(h_fr, eb)


h2d = ROOT.TH2F('photon_fakes', '', 54, 30, 300, 2, array('d', [0., 1.4442, 2.5]))

for i in xrange(1, h2d.GetNbinsX() + 1):
    x = h2d.GetXaxis().GetBinCenter(i)
    h2d.SetBinContent(i, 1, res['barrel_%s' % args.channel].GetBinContent(res['barrel_%s' % args.channel].GetXaxis().FindFixBin(x)))
    h2d.SetBinContent(i, 2, res['endcap_%s' % args.channel].GetBinContent(res['endcap_%s' % args.channel].GetXaxis().FindFixBin(x)))

ROOT.gDirectory.WriteTObject(h2d, 'photon_fakes')

fout.Close()