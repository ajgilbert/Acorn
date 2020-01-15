import ROOT
import argparse
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(True)
plot.ModTDRStyle()

parser = argparse.ArgumentParser()
parser.add_argument('-i', help="input.root:tree")
parser.add_argument('-d', default='', help='Draw expression')
parser.add_argument('-b', default='', help='Binning')
parser.add_argument('-s', default='', help='Selection/weights')
parser.add_argument('-o', default='', help='output.root:hist')
parser.add_argument('-x', default='', help='x-axis title')
args = parser.parse_args()
# file:tree:var:binning:sel
fin = ROOT.TFile(args.i.split(':')[0])
tree = fin.Get(args.i.split(':')[1])

tree.Draw('%s>>h%s' % (args.d, args.b), args.s)

h = ROOT.gDirectory.Get('h')
h.SetName(args.o.split(':')[1])

canv = ROOT.TCanvas(args.o.split(':')[1], args.o.split(':')[1])
pads = plot.OnePad()
h.Draw('HIST')
h.GetXaxis().SetTitle(args.x)
canv.Print('.pdf')
canv.Print('.png')

fout = ROOT.TFile(args.o.split(':')[0], 'RECREATE')
h.Write()
fout.Close()

fin.Close()
