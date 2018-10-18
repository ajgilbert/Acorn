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

fin = ROOT.TFile('output_2016_photon_fakes.root')
hists = TDirToNode(fin)

var = 'p0_pt'
# var = 'p0_chiso'

def RebinHist(hist, bins):
    x = hist.Rebin(len(bins) - 1, "", array('d', bins))
    x.Copy(hist)


fout = ROOT.TFile('output_2016_photon_fakes_ratios.root', 'RECREATE')

res = {}

for eb in ['barrel', 'endcap']:
    for sel in ['%s_iso_t_sig_l' % eb, '%s_iso_l_sig_l' % eb, '%s_iso_l_sig_t' % eb]:
        node = hists[sel][var]
        if eb == 'barrel':
            newbins = [30, 35, 40, 50, 60, 80, 100, 300]
        else:
            newbins = [30, 40, 60, 300]

        # newbins = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        node.ForEach(lambda x: RebinHist(x, newbins))
        node['data_sub'] = node['data_obs'] - node['Total_R']
    h_fr = hists['%s_iso_l_sig_t' % eb][var]['data_sub'].Clone()
    h_fr.Divide( hists['%s_iso_l_sig_l' % eb][var]['data_sub'])

    h_fr_mc = hists['%s_iso_l_sig_t' % eb][var]['W_F'].Clone()
    h_fr_mc.Divide( hists['%s_iso_l_sig_l' % eb][var]['W_F'])

    plot.Set(h_fr, MarkerSize=0.5)
    canv = ROOT.TCanvas('photon_fakes_%s' % eb, 'photon_fakes_%s' % eb)
    pads = plot.OnePad()
    h_fr.Draw('E')
    plot.Set(h_fr_mc, LineColor=2, MarkerColor=2)
    # h_fr_mc.Draw('SAMEE')
    h_fr.SetMaximum(2.0)
    h_fr.SetMinimum(0.0)
    canv.Print('.pdf')
    canv.Print('.png')
    fout.cd()
    res[eb] = h_fr
    ROOT.gDirectory.WriteTObject(h_fr, eb)


h2d = ROOT.TH2F('photon_fakes', '', 54, 30, 300, 2, array('d', [0., 1.4442, 2.5]))

for i in xrange(1, h2d.GetNbinsX() + 1):
    x = h2d.GetXaxis().GetBinCenter(i)
    h2d.SetBinContent(i, 1, res['barrel'].GetBinContent(res['barrel'].GetXaxis().FindFixBin(x)))
    h2d.SetBinContent(i, 2, res['endcap'].GetBinContent(res['endcap'].GetXaxis().FindFixBin(x)))

ROOT.gDirectory.WriteTObject(h2d, 'photon_fakes')

fout.Close()