import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
import ROOT
import argparse
import json
import os
import fnmatch
from copy import deepcopy
from array import array

ROOT.TH1.SetDefaultSumw2(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()


parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='Output of PostFitShapes or PostFitShapesFromWorkspace, specified as FILE:BIN')
# parser.add_argument('--output', '-o', default='.', help='top level output directory')
# parser.add_argument('--channel', '-c', default='mt', choices=['mt', 'et', 'em', 'tt', 'mm', 'generic', 'wgamma'], help='Channel')
# parser.add_argument('--x-title', default='', help='x-axis variable, without GeV')
# parser.add_argument('--logy', action='store_true')
# parser.add_argument('--y-min', type=float, default=1)
# parser.add_argument('--title', default='')
parser.add_argument('--layout-file', '-l', default='layouts.json')
parser.add_argument('--cut-gt', default=None)
parser.add_argument('--cut-lt', default=None)
parser.add_argument('--overflow', action='store_true')
parser.add_argument('--layout', default='data_fakes')

args = parser.parse_args()

with open(args.layout_file) as jsonfile:
    layouts = json.load(jsonfile)

layout = layouts[args.layout]

filename, dirfilter = args.input.split(':')
print filename
file = ROOT.TFile(filename)

node = TDirToNode(file)

made_dirs = set()

for path, subnode in node.ListNodes(withObjects=True):
    if not fnmatch.fnmatch(path, dirfilter):
        continue
    # for now work on the assumption that the last component of the path will be the actual filename
    print path
    split_path = path.split('/')[:-1]
    name = path.split('/')[-1]
    # target_dir = os.path.join(args.output, *split_path)
    hists = {}
    for opath, objname, obj in subnode.ListObjects(depth=0):
        hists[objname] = obj

    for ele in layout:
        htemp = hists[ele['entries'][0]].Clone()
        for oth in ele['entries'][1:]:
            htemp.Add(hists[oth])
        bin_lo = 1
        bin_hi = htemp.GetNbinsX()
        if args.overflow:
            bin_lo = 0
            bin_hi += 1
        if args.cut_gt != None:
            bin_lo = htemp.GetXaxis().FindFixBin(float(args.cut_gt))
            # print htemp.GetXaxis().GetBinLowEdge(bin_lo)
        if args.cut_lt != None:
            bin_hi = htemp.GetXaxis().FindFixBin(float(args.cut_lt))
            # print htemp.GetXaxis().GetBinUpEdge(bin_hi)
        integral = htemp.Integral(bin_lo, bin_hi)
        print '%-30s : %20.1f' % (ele['name'], integral)
