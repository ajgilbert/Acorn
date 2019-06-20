import ROOT
import math
import argparse
import sys
import os
from array import array

#ROOT.gROOT.SetBatch(ROOT.kTrue)
#ROOT.TH1.AddDirectory(False)

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help="Input file")
parser.add_argument('--outname','-o', help="Name for output file")
parser.add_argument('--pull','-p', action='store_true', help="Plot pull", default=False)
parser.add_argument('-r', help="injected signal strength",type=float) 
parser.add_argument('--custom_x_range', help='Fix x axis range', action='store_true', default=False)
parser.add_argument('--xmin',default="-0.1",help="x min")
parser.add_argument('--xmax',default="0.1",help="x max")


args = parser.parse_args()
xmin = -0.05
xmax = 0.05
if args.pull:
  xmin = -4.0
  xmax = 4.0
if args.custom_x_range:
  xmin = float(args.xmin)
  xmax = float(args.xmax)
    

inf = ROOT.TFile.Open(args.file)
plot_tree = inf.Get("tree_fit_sb")
ROOT.gStyle.SetOptStat(1110)

outputhisto = ROOT.TH1F("outputhisto","outputhisto",50,xmin,xmax)

outputhisto.GetYaxis().SetTitle("nEntries")
outputhisto.GetXaxis().SetTitle("Abs bias")
if args.pull:
  outputhisto.GetXaxis().SetTitle("Pull")


c1 = ROOT.TCanvas()
c1.cd()

if args.pull:
  plot_tree.Draw("(r-%f)/rErr>>outputhisto"%float(args.r),"fit_status>-0.5")
else:
  plot_tree.Draw("(r-%f)>>outputhisto"%float(args.r),"fit_status>-0.5")

outputhisto.SetTitle("")
addout = "%s"%(args.r)

outname = args.outname+"_"+addout

if args.pull:
  outname = outname+"_pull"

outputhisto.Fit("gaus")

c1.SaveAs("%(outname)s.pdf"%vars())
