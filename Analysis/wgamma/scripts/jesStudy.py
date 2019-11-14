#!/usr/bin/env python
import ROOT
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.plottemplates import *
import argparse
import sys
import math

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.SetDefaultSumw2(True)
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


def PlotCorrMatrix(N, corr, jecs):
    plot.ModTDRStyle(width=750, height=600, l = 0.20, b = 0.20, r = 0.23, t=0.08)
    plot.SetCorrMatrixPalette()
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetMarkerSize(1.0)
    ROOT.gStyle.SetPaintTextFormat('.2f')
    hist = ROOT.TH2F('cor', '', N, 0, N, N, 0, N)
    for i, ip in enumerate(jecs):
        hist.GetXaxis().SetBinLabel(i + 1, ip)
        hist.GetYaxis().SetBinLabel(N - i, ip)
        for j, jp in enumerate(jecs):
            hist.SetBinContent(i + 1, N - j, corr[i][j])
    hist.GetXaxis().LabelsOption('v')
    canv = ROOT.TCanvas('jes_corr', 'jes_corr')
    pads = plot.OnePad()
    ROOT.gStyle.SetTextFont(42)
    # hist.SetMarkerSize(args.marker_size)
    hist.Draw('COLZ')
    canv.Print('.pdf')
    canv.Print('.png')


def SimpleUnroll(hist2d):
    nbins = hist2d.GetNbinsX() * hist2d.GetNbinsY()
    hist = ROOT.TH1F(hist2d.GetName(), '', nbins, 0, nbins)
    ib = 1
    for iy in xrange(1, hist2d.GetNbinsY() + 1):
        for ix in xrange(1, hist2d.GetNbinsX() + 1):
            hist.SetBinContent(ib, hist2d.GetBinContent(ix, iy))
            hist.SetBinError(ib, hist2d.GetBinError(ix, iy))
            hist.GetXaxis().SetBinLabel(ib, '[%.0f,%.0f]' % (hist2d.GetXaxis().GetBinLowEdge(ix), hist2d.GetXaxis().GetBinUpEdge(ix)))
            ib += 1
    hist.GetXaxis().LabelsOption('v')
    return hist

def DrawUnrolled(h_nominal, jecs, hists, label):
    plot.ModTDRStyle(width=1000, b=0.2)
    newhists = {
        'nominal': SimpleUnroll(h_nominal)
    }
    for jec in jecs:
        newhists[jec] = SimpleUnroll(hists[jec])
    plotcfg = DEFAULT_CFG

    # for h in newhists:
    #     newhists[h].Scale(1.0 / newhists[h].Integral())
    layout=[
        {'name': 'nominal', 'legend': 'Nominal', 'marker_size': 0}
    ]
    ratios = [
    ]
    for jec in jecs:
        layout.append({'name': jec, 'legend': jec, 'marker_size': 0})
        ratios.append({'num': jec, 'den': 'nominal', 'type': 'noerror'})

    plotcfg.update({
        'type': 'multihist',
        'legend_pos': [0.55, 0.8, 0.90, 0.91],
        'ratio_pad_frac': 0.45,
        'logy': True,
        'logy_min': 1E-4,
        'global_hist_opts': {
            'draw_opts': 'L',
            'legend_opts': 'L',
            'marker_size': 0.0,
            'line_width': 1
        },
        })
    MakeMultiHistPlot(label,
                      outdir='.',
                      hists=newhists,
                      cfg=UpdateDict(plotcfg, {
                          'x_title': 'p_{T},#eta bin',
                          'ratio_y_title': 'Shifted/Nominal',
                          'ratio_y_range': [0.95, 1.05],
                          }),
                      layout=layout,
                      ratios=ratios
    )




def FindLargestCorrelation(N, corr, labels, veto=[]):
    largest = 0.0
    l1 = 0
    l2 = 0
    for i in xrange(N):
        for j in xrange(i + 1, N):
            if labels[i] in veto or labels[j] in veto:
                continue
            if abs(corr[i][j]) > largest:
                largest = abs(corr[i][j])
                l1 = i
                l2 = j
    print '>> Largest correlation between %s (%i) and %s (%i): %f' % (labels[l1], l1, labels[l2], l2, largest)
    return [labels[l1], labels[l2]]


def FindOthers(N, corr, labels, selected, threshold=0.9):
    res = list()
    for i in xrange(N):
        if labels[i] in selected:
            print '>> Skip %s, already selected' % labels[i]
            continue
        ok = True
        for sel in selected:
            idx = labels.index(sel)
            c = corr[i][idx]
            if abs(c) < threshold:
                ok = False
        if ok:
            res.append(labels[i])
    return res

# parser = argparse.ArgumentParser()
# args = parser.parse_args()

plot.ModTDRStyle(width=750, height=600, l=0.17, b=0.15, r=0.23, t=0.08)

# plot.SetCorrMatrixPalette()
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(1.0)
ROOT.gStyle.SetPaintTextFormat('.2f')

fin = ROOT.TFile(sys.argv[1])
subdir = sys.argv[2]

hists = {}
hists['nominal'] = fin.Get('%s/nominal' % subdir)
for jec in jecs:
    hists[jec] = fin.Get('%s/%s' % (subdir, jec))

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
        if TH2Mag(hists_normed[jecs[i]]) <= 0.:
            print jecs[i]
        if TH2Mag(hists_normed[jecs[j]]) <= 0.:
            print jecs[j]
        corr[i][j] = x[i][j] / (TH2Mag(hists_normed[jecs[i]]) * TH2Mag(hists_normed[jecs[j]]))
# x.Print()
corr.Print()

correlations = []
for i in xrange(N):
    for j in xrange(i + 1, N):
        if abs(corr[i][j]) > 0.8:
            correlations.append((jecs[i], jecs[j], corr[i][j]))

correlations.sort(key=lambda x: abs(x[2]), reverse=True)
# print correlations
for i in xrange(len(correlations)):
    print '%i  %30s (%i) %30s (%i) %+.3f' % (i, correlations[i][0], jecs.index(correlations[i][0]), correlations[i][1], jecs.index(correlations[i][1]), correlations[i][2])



PlotCorrMatrix(N, corr, jecs)

testgroups = [
    [
        'AbsoluteStat',
        'AbsoluteScale',
        'AbsoluteMPFBias',
        'RelativeStatFSR'
    ],
    [
        'Fragmentation',
        'SinglePionECAL',
        'SinglePionHCAL',
        'FlavorQCD',
        'RelativeFSR',
        'PileUpDataMC',
    ],
    [
        'TimePtEta',
        'RelativeJEREC2',
        'RelativePtEC2',
        'RelativeStatEC'
    ],
    [
        'RelativeJERHF',
        'RelativePtHF',
        'RelativeStatHF',
    ],
    [
        'RelativeJEREC1',
        'RelativePtEC1',
        'PileUpPtEC1'
    ],
    [
        'PileUpPtRef',
        'PileUpPtHF'
    ],
    [
        'AbsoluteSample',
        'RelativeSample'
    ]
]

groups_htt = [
    ["SinglePionECAL", "SinglePionHCAL", "AbsoluteMPFBias", "AbsoluteScale", "AbsoluteStat", "Fragmentation", "FlavorQCD", "TimePtEta", "PileUpDataMC", "RelativeFSR", "RelativeStatFSR", "PileUpPtRef"],
    ["PileUpPtEC1", "PileUpPtEC2", "PileUpPtBB", "RelativeJEREC1", "RelativeJEREC2", "RelativePtEC1", "RelativePtEC2", "RelativeStatEC", "RelativePtBB"],
    ["RelativeStatHF", "RelativePtHF", "PileUpPtHF", "RelativeJERHF"],
    ["RelativeBal"],
    ["RelativeSample"],
]

for igrp, grp in enumerate(groups_htt):
    # DrawUnrolled(hists['nominal'], grp, hists, 'group%i' % igrp)
    print 'Testing group: %s' % grp
    for i in range(len(grp)):
        for j in range(i + 1, len(grp)):
            this_corr = corr[jecs.index(grp[i])][jecs.index(grp[j])]
            if abs(this_corr) < 0.8:
                print '>> %s-%s = %f' % (grp[i], grp[j], this_corr)



# best_list = FindLargestCorrelation(N, corr, jecs)
# add_list = FindOthers(N, corr, jecs, best_list, 0.9)
# print add_list

# best_list = FindLargestCorrelation(N, corr, jecs, veto=(best_list + add_list))

# add_list = FindOthers(N, corr, jecs, best_list, 0.8)
# print add_list

# eigen = ROOT.TMatrixDSymEigen(x)
# evals = eigen.GetEigenValues()
# evecs = eigen.GetEigenVectors()

# evals.Print()
# evecs.Print()



# print hists
