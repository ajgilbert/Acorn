import ROOT
import argparse
import json
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
from Acorn.Analysis.plottemplates import *


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)

parser = argparse.ArgumentParser()
# parser.add_argument('--selection', default='eft_region')
parser.add_argument('--output', default='diff_xsec_breakdown')
parser.add_argument('--scheme', default='phi_f_binned')
# parser.add_argument('--charge', default='x')
# parser.add_argument('--split-unc', default=0, type=int)
args = parser.parse_args()

chg = 'x'
chg_labels = {
    'p': '+',
    'n': '-',
    'x': '#pm'
}

plot.ModTDRStyle()

canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()

legend = ROOT.TLegend(0.55, 0.8, 0.95, 0.92, '', 'NBNDC')
legend.SetNColumns(2)

xsec_results = {}
with open('%s.json' % args.scheme) as jsonfile:
    xsec_results = json.load(jsonfile)['xsec2D']

bin_edges = array('d', [30,50,70,100,150,200,300,500,800,1200])

sources = ['Error', 'Stat', 'Lumi', 'Expt', 'MCStats', 'Th']
labels = ['Total', 'Stat', 'Lumi', 'Expt', 'MCStats', 'Theory']
cols = [1, 2, 4, 6, 28, 8, 9, 13]
hists = []
for i, s in enumerate(sources):
    hist = ROOT.TH1D(s, s, len(bin_edges) - 1, bin_edges)
    for ib in xrange(hist.GetNbinsX()):
        av_err = 0.5 * (xsec_results['r_x_%i' % ib]['%sHi' % sources[i]] - xsec_results['r_x_%i' % ib]['%sLo' % sources[i]])
        rel_av_err = av_err / abs(xsec_results['r_x_%i' % ib]['Val'])
        hist.SetBinContent(ib + 1, rel_av_err)
    plot.Set(hist, LineColor=cols[i], LineWidth=3)
    if i == 0:
        plot.Set(hist, LineColor=cols[i], LineWidth=4)
    legend.AddEntry(hist, labels[i], 'L')
    hists.append(hist)
    hist.Draw('HISTSAME')

axis = hists[0]
axis.GetXaxis().SetTitle('Photon p_{T} (GeV)')
axis.GetYaxis().SetTitle('Rel. uncertainty on #sigma_{i}')

axis.SetMinimum(0.01)
pads[0].SetLogy()
pads[0].SetGrid(0, 1)
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.25)


legend.Draw()

pads[0].cd()
plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)

plot.DrawTitle(pads[-1], '136.9 fb^{-1} (13 TeV)', 3)
plot.DrawTitle(pads[0], 'W^{%s}(l^{%s}#nu)#gamma' % (chg_labels[chg], chg_labels[chg]), 1)

canv.Print('.png')
canv.Print('.pdf')
