import ROOT
import argparse
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(True)
plot.ModTDRStyle()

parser = argparse.ArgumentParser()
parser.add_argument('--get-quantiles', default=None, help="input.root:hist")
parser.add_argument('--apply-quantiles', default=None, help="input.root:hist")
# parser.add_argument('-d', default='', help='Draw expression')
# parser.add_argument('-b', default='', help='Binning')
# parser.add_argument('-s', default='', help='Selection/weights')
# parser.add_argument('-o', default='', help='output.root:hist')
# parser.add_argument('-x', default='', help='x-axis title')
args = parser.parse_args()

filename = args.get_quantiles
# file:tree:var:binning:sel
fin = ROOT.TFile(filename.split(':')[0])
h = fin.Get(filename.split(':')[1])

integral = h.Integral()

edges = []
quantiles = []
fracs = []
running_sum = 0.

# h.Print('all')

for ix in xrange(1, h.GetNbinsX() + 1):
    xlow = h.GetXaxis().GetBinLowEdge(ix)
    edges.append(xlow)
    quantiles.append(running_sum / integral)
    fracs.append(h.GetBinContent(ix) / integral)
    running_sum += h.GetBinContent(ix)
    if ix == h.GetNbinsX():
        edges.append(h.GetXaxis().GetBinUpEdge(ix))
        quantiles.append(running_sum / integral)
        fracs.append(h.GetBinContent(ix) / integral)

for X in zip(edges, quantiles, fracs):
    print X

fin.Close()

res_q = array('d', [0.] * len(quantiles))

if args.apply_quantiles is not None:
    filename = args.apply_quantiles
    fin = ROOT.TFile(filename.split(':')[0])
    h = fin.Get(filename.split(':')[1])

    h.GetQuantiles(len(quantiles), res_q, array('d', quantiles))

    for X in zip(res_q, quantiles):
        print X
