#!/usr/bin/env python
import ROOT
import CombineHarvester.CombineTools.plotting as plot
import argparse
import sys
import math

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

jecs = [
    "AbsoluteStat",
    "AbsoluteScale",
    "AbsoluteMPFBias",
    "AbsoluteSample",
    "Fragmentation",
    "SinglePionECAL",
    "SinglePionHCAL",
    "FlavorQCD",
    "TimePtEta",
    "RelativeJEREC1",
    "RelativeJEREC2",
    "RelativeJERHF",
    "RelativePtBB",
    "RelativePtEC1",
    "RelativePtEC2",
    "RelativePtHF",
    "RelativeBal",
    "RelativeSample",
    "RelativeFSR",
    "RelativeStatFSR",
    "RelativeStatEC",
    "RelativeStatHF",
    "PileUpDataMC",
    "PileUpPtRef",
    "PileUpPtBB",
    "PileUpPtEC1",
    "PileUpPtEC2",
    "PileUpPtHF"
]


def Merge(jecs, jec1, jec2):
    newjecs = []
    for jec in jecs:
        if jec == jec1:
            newjecs.append('%s+%s' % (jec1, jec2))
        elif jec == jec2:
            continue
        else:
            newjecs.append(jec)
    return newjecs



def Add(hist2d, val):
    for ix in xrange(1, hist2d.GetNbinsX() + 1):
        for iy in xrange(1, hist2d.GetNbinsY() + 1):
            hist2d.SetBinContent(ix, iy, hist2d.GetBinContent(ix, iy) + val)


def TH2DotProd(h1, h2):
    val = 0.
    for ix in xrange(1, h1.GetNbinsX() + 1):
        for iy in xrange(1, h1.GetNbinsY() + 1):
            val += h1.GetBinContent(ix, iy) * h2.GetBinContent(ix, iy)
    return val


def TH2Mag(h1):
    val = 0.
    for ix in xrange(1, h1.GetNbinsX() + 1):
        for iy in xrange(1, h1.GetNbinsY() + 1):
            val += h1.GetBinContent(ix, iy) * h1.GetBinContent(ix, iy)
    return math.sqrt(val)

# parser = argparse.ArgumentParser()
# args = parser.parse_args()

plot.ModTDRStyle(width=750, height=600, l=0.17, b=0.15, r=0.23, t=0.08)

# plot.SetCorrMatrixPalette()
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(1.0)
ROOT.gStyle.SetPaintTextFormat('.2f')

fin = ROOT.TFile(sys.argv[1])

hists = {}
hists['nominal'] = fin.Get('nominal')
for jec in jecs:
    hists[jec] = fin.Get(jec)

hists_normed = {}
for i, jec in enumerate(jecs):
    hists_normed[jec] = hists[jec].Clone()
    hists_normed[jec].Divide(hists['nominal'])
    Add(hists_normed[jec], -1.0)
    print i, jec
    # hists_normed[jec].Print('range')

N = len(jecs)

x = ROOT.TMatrixDSym(N)
corr = ROOT.TMatrixDSym(N)

for i in xrange(N):
    for j in xrange(N):
        x[i][j] = TH2DotProd(hists_normed[jecs[i]], hists_normed[jecs[j]])
        corr[i][j] = x[i][j] / (TH2Mag(hists_normed[jecs[i]]) * TH2Mag(hists_normed[jecs[j]]))
# x.Print()
corr.Print()

eigen = ROOT.TMatrixDSymEigen(x)
evals = eigen.GetEigenValues()
evecs = eigen.GetEigenVectors()

evals.Print()
evecs.Print()



# print hists
