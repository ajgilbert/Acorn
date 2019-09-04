import ROOT
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
import argparse
import json

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

parser = argparse.ArgumentParser()
parser.add_argument(
    'input', nargs='+', help="""Input json files""")
parser.add_argument(
    '--show', default='1sig,2sig')
parser.add_argument(
    '--output', '-o', default='limit', help="""Name of the output
    plot without file extension""")
parser.add_argument(
    '--x-title', default='m_{H} (GeV)', help="""Title for the x-axis""")
parser.add_argument(
    '--y-title', default=None, help="""Title for the y-axis""")
parser.add_argument(
    '--y-range', default='-10,10', help="""Range for the y-axis""")
parser.add_argument('--pt-bins', default=None)

args = parser.parse_args()

show = args.show.split(',')

def DrawAxisHists(pads, axis_hists, def_pad=None, pt_bins=None):
    for i, pad in enumerate(pads):
        pad.cd()
        if pt_bins is not None:
            axis_hists[i].GetXaxis().SetLimits(pt_bins[0], pt_bins[-1])
        axis_hists[i].Draw('AXIS')
        axis_hists[i].Draw('AXIGSAME')
    if def_pad is not None:
        def_pad.cd()

pt_bins = [150, 210, 300, 420, 600, 850, 1200]

bin_centres = []
bin_e_lo = []
bin_e_hi = []

for i in xrange(len(pt_bins) - 1):
    bin_centres.append((pt_bins[i+1] + pt_bins[i]) / 2.)
    bin_e_hi.append((pt_bins[i + 1] - pt_bins[i]) / 2.)
    bin_e_lo.append((pt_bins[i + 1] - pt_bins[i]) / 2.)

canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()

for padx in pads:
    # Use tick marks on oppsite axis edges
    plot.Set(padx, Tickx=1, Ticky=1, Logx=False)

graphs = []
graph_sets = []

legend = plot.PositionedLegend(0.30, 0.20, 3, 0.015)
# plot.Set(legend, NColumns=2)

axis = None
altStyle = True

# Process each input argument
for src in args.input:
    splitsrc = src.split(':')
    file = splitsrc[0]

    with open(file) as jsonfile:
        js = json.load(jsonfile)

    m2 = []
    m1 = []
    y = []
    p1 = []
    p2 = []
    print js['c3w_bin_%i' % i]['c3w']

    for i in xrange(len(bin_centres)):
        y.append(js['c3w_bin_%i' % i]['c3w']['Val'])
        if altStyle:
            m2.append(js['c3w_bin_%i' % i]['c3w']['2sig_ErrorLo'] + y[i])
            m1.append(js['c3w_bin_%i' % i]['c3w']['ErrorLo'] + y[i])
            p1.append(js['c3w_bin_%i' % i]['c3w']['ErrorHi'] + y[i])
            p2.append(js['c3w_bin_%i' % i]['c3w']['2sig_ErrorHi'] + y[i])
        else:
            m2.append(js['c3w_bin_%i' % i]['c3w']['2sig_ErrorLo'] * -1.)
            m1.append(js['c3w_bin_%i' % i]['c3w']['ErrorLo'] * -1.)
            p1.append(js['c3w_bin_%i' % i]['c3w']['ErrorHi'])
            p2.append(js['c3w_bin_%i' % i]['c3w']['2sig_ErrorHi'])

    graph_set = {}

    if altStyle:
        graph_set['nominal'] = ROOT.TGraph(len(bin_centres), array('d', pt_bins[1:]), array('d', y))
        graph_set['1sig_p'] = ROOT.TGraph(len(bin_centres), array('d', pt_bins[1:]), array('d', p1))
        graph_set['1sig_m'] = ROOT.TGraph(len(bin_centres), array('d', pt_bins[1:]), array('d', m1))
        graph_set['2sig_p'] = ROOT.TGraph(len(bin_centres), array('d', pt_bins[1:]), array('d', p2))
        graph_set['2sig_m'] = ROOT.TGraph(len(bin_centres), array('d', pt_bins[1:]), array('d', m2))
        col = int(splitsrc[1])
        plot.Set(graph_set['nominal'], MarkerColor=col, LineColor=col, LineWidth=2, MarkerSize=0.8)
        plot.Set(graph_set['1sig_p'], MarkerColor=col, LineColor=col, LineWidth=2, MarkerSize=0.5, LineStyle=2)
        plot.Set(graph_set['1sig_m'], MarkerColor=col, LineColor=col, LineWidth=2, MarkerSize=0.5, LineStyle=2)
        plot.Set(graph_set['2sig_p'], MarkerColor=col, LineColor=col, LineWidth=2, MarkerSize=0.5, LineStyle=9)
        plot.Set(graph_set['2sig_m'], MarkerColor=col, LineColor=col, LineWidth=2, MarkerSize=0.5, LineStyle=9)
    else:
        graph_set['nominal'] = ROOT.TGraphAsymmErrors(len(bin_centres), array('d', pt_bins[1:]), array('d', y), array('d', [0]*len(bin_centres)), array('d', [0]*len(bin_centres)), array('d', [0]*len(bin_centres)), array('d', [0]*len(bin_centres)))
        graph_set['1sig'] = ROOT.TGraphAsymmErrors(len(bin_centres), array('d', pt_bins[1:]), array('d', y), array('d', [0]*len(bin_centres)), array('d', [0]*len(bin_centres)), array('d', m1), array('d', p1))
        graph_set['2sig'] = ROOT.TGraphAsymmErrors(len(bin_centres), array('d', pt_bins[1:]), array('d', y), array('d', [0]*len(bin_centres)), array('d', [0]*len(bin_centres)), array('d', m2), array('d', p2))
        plot.Set(graph_set['2sig'], LineColor=2, LineWidth=5, MarkerSize=0)
        plot.Set(graph_set['1sig'], LineColor=4, LineWidth=10, MarkerSize=0)

    graph_sets.append(graph_set)

    if axis is None:
        axis = plot.CreateAxisHists(len(pads), graph_sets[-1].values()[0], True)
        DrawAxisHists(pads, axis, pads[0], pt_bins)
        axis[0].GetXaxis().SetLimits(pt_bins[1] - (0.05 * (pt_bins[-1] - pt_bins[1])), pt_bins[-1] + (0.05 * (pt_bins[-1] - pt_bins[1])))
        axis[0].SetMinimum(float(args.y_range.split(',')[0]))
        axis[0].SetMaximum(float(args.y_range.split(',')[1]))
        line = ROOT.TLine()
        plot.DrawHorizontalLine(pads[0], line, 0.)

    if altStyle:
        graph_set['nominal'].Draw('PLSAME')
        legend.AddEntry(graph_set['nominal'], '%s' % splitsrc[2], 'LP')
        if '1sig' in args.show:
            graph_set['1sig_p'].Draw('PLSAME')
            graph_set['1sig_m'].Draw('PLSAME')
            legend.AddEntry(graph_set['1sig_p'], '#pm 1#sigma %s' % splitsrc[2], 'L')
        if '2sig' in args.show:
            graph_set['2sig_p'].Draw('PLSAME')
            graph_set['2sig_m'].Draw('PLSAME')
            legend.AddEntry(graph_set['2sig_p'], '#pm 2#sigma %s' % splitsrc[2], 'L')
    else:
        graph_set['2sig'].Draw('PZSAME')
        graph_set['1sig'].Draw('PZSAME')
        graph_set['nominal'].Draw('PZSAME')



axis[0].GetYaxis().SetTitle('C_{3W} [TeV^{-2}]')
if args.y_title is not None:
    axis[0].GetYaxis().SetTitle(args.y_title)
axis[0].GetXaxis().SetTitle(args.x_title)
axis[0].GetXaxis().SetLabelOffset(axis[0].GetXaxis().GetLabelOffset()*2)

# pads[0].SetLogy(1)

# plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)
legend.Draw()

pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

# # CMS logo
plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)
# # plot.DrawTitle(pads[0], '0.1 fb^{-1} (13 TeV)', 3)

canv.Print('.png')
canv.Print('.pdf')
