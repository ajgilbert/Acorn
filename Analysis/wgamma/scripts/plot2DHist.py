#!/usr/bin/env python
import ROOT
import CombineHarvester.CombineTools.plotting as plot
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def ColorLabels(text):
    return text


parser = argparse.ArgumentParser()
parser.add_argument('--output', '-o', help='output name')
parser.add_argument('--input', '-i', help='Main input file')
parser.add_argument('--label', '-l', help='Main input file')
parser.add_argument('--label-size', type=float, default=0.04)
parser.add_argument('--marker-size', type=float, default=1.5)
parser.add_argument('--print-twiki', action='store_true')
parser.add_argument('--titles', default='x,y,z')
parser.add_argument(
    '--subline', default='Internal', help='text to add next to cms logo')

args = parser.parse_args()
plot.ModTDRStyle(width=750, height=600, l = 0.17, b = 0.15, r = 0.23, t=0.08)

# plot.SetCorrMatrixPalette()
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(1.0)
ROOT.gStyle.SetPaintTextFormat('.2f')

x_title, y_title, z_title = args.titles.split(',')

inputs = args.input.split(':')
fin = ROOT.TFile(inputs[0])
hist = fin.Get(inputs[1])

hist.GetXaxis().SetTitle(x_title)
hist.GetYaxis().SetTitle(y_title)
hist.GetZaxis().SetTitle(z_title)
hist.GetXaxis().SetLabelFont(42)
hist.GetYaxis().SetLabelFont(42)
hist.GetXaxis().SetTickLength(0)
hist.GetYaxis().SetTickLength(0)
hist.GetXaxis().SetLabelSize(args.label_size)
hist.GetYaxis().SetLabelSize(args.label_size)
hist.GetZaxis().SetLabelSize(0.045)
# hist.GetZaxis().SetTitle('#rho')
hist.SetContour(255)

canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()


ROOT.gStyle.SetTextFont(42)
hist.SetMarkerSize(args.marker_size)
hist.Draw('COLZTEXT')


plot.DrawCMSLogo(pads[0], 'CMS', args.subline, 0, 0.17, 0.035, 1.2)

if args.label is not None:
    plot.DrawTitle(pads[0], args.label, 3)


canv.Print('.pdf')
canv.Print('.png')
