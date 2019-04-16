import ROOT
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

parser = argparse.ArgumentParser()
parser.add_argument(
    'input', nargs='+', help="""Input json files""")
parser.add_argument(
    '--show', default='exp,obs')
parser.add_argument(
    '--output', '-o', default='limit', help="""Name of the output
    plot without file extension""")
parser.add_argument(
    '--x-title', default='m_{H} (GeV)', help="""Title for the x-axis""")
parser.add_argument(
    '--y-title', default=None, help="""Title for the y-axis""")
parser.add_argument(
    '--limit-on', default='#sigma/#sigma_{SM}', help="""Shortcut for setting the y-axis label""")
parser.add_argument('--pt-bins', default=None)

args = parser.parse_args()


def DrawAxisHists(pads, axis_hists, def_pad=None, pt_bins=None):
    for i, pad in enumerate(pads):
        pad.cd()
        if pt_bins is not None:
            axis_hists[i].GetXaxis().SetLimits(pt_bins[0], pt_bins[-1])
        axis_hists[i].Draw('AXIS')
        axis_hists[i].Draw('AXIGSAME')
    if def_pad is not None:
        def_pad.cd()


def ConvertToHist(gr, bins):
    N = gr.GetN()
    hist = ROOT.TH1F('', '', N, array('d', pt_bins[0:N + 1]))
    for i in xrange(N):
        hist.SetBinContent(i + 1, gr.GetY()[i])
    gr.Print()
    hist.Print("range")
    return hist


def AdjustGraph(gr, bins):
    if not isinstance(gr, ROOT.TGraphAsymmErrors):
        new_gr = ROOT.TGraphAsymmErrors(gr.GetN(), gr.GetX(), gr.GetY())
    else:
        new_gr = gr.Clone()
    # plot.Set(new_gr, LineColor=gr.GetLineColor(), LineWidth=gr.GetLineWidth(), LineStyle=gr.GetLineStyle(),
    #                  MarkerSize=gr.GetMarkerSize(), MarkerColor=gr.GetMarkerColor(), MarkerStyle=gr.GetMarkerStyle(),
    #                  FillColor=gr.GetFillColor(), FillStyle=gr.GetFillStyle())
    for i in range(gr.GetN()):
        new_gr.GetX()[i] = (bins[i+1] + bins[i]) / 2.
        new_gr.GetEXhigh()[i] = (bins[i + 1] - bins[i]) / 2.
        new_gr.GetEXlow()[i] = (bins[i + 1] - bins[i]) / 2.
    return new_gr


pt_bins = [150, 210, 300, 420, 600, 850, 1200]

canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()

for padx in pads:
    # Use tick marks on oppsite axis edges
    plot.Set(padx, Tickx=1, Ticky=1, Logx=False)

graphs = []
graph_sets = []

legend = plot.PositionedLegend(0.45, 0.10, 3, 0.015)
plot.Set(legend, NColumns=2)

axis = None

# Process each input argument
for src in args.input:
    splitsrc = src.split(':')
    file = splitsrc[0]
    # limit.json => Draw as full obs + exp limit band
    if len(splitsrc) == 1:
        graph_sets.append(plot.StandardLimitsFromJSONFile(file, args.show.split(',')))
        for gr_name in graph_sets[-1]:
            graph_sets[-1][gr_name] = AdjustGraph(graph_sets[-1][gr_name], pt_bins)
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), graph_sets[-1].values()[0], True)
            DrawAxisHists(pads, axis, pads[0], pt_bins)
        plot.StyleLimitBand(graph_sets[-1], overwrite_style_dict={
            'exp0': {'MarkerSize': 0}
            })
        plot.DrawLimitBand(pads[0], graph_sets[-1], legend=legend, legend_overwrite={
            'exp2': {'DrawStyle': '2SAME'},
            'exp1': {'DrawStyle': '2SAME'},
            'obs': {'DrawStyle': 'PZSAME'},
            'exp0': {'DrawStyle': 'PZSAME'}
            })
        pads[0].RedrawAxis()
        pads[0].RedrawAxis('g')
        pads[0].GetFrame().Draw()


    # limit.json:X => Draw a single graph for entry X in the json file
    # 'limit.json:X:Title="Blah",LineColor=4,...' =>
    # as before but also apply style options to TGraph
    elif len(splitsrc) >= 2:
        settings = {}
        settings['Title'] = src
        graphs.append(plot.LimitTGraphFromJSONFile(file, splitsrc[1]))
        graphs[-1] = AdjustGraph(graphs[-1], pt_bins)
        if len(splitsrc) >= 3:
            settings.update({x.split('=')[0]: eval(x.split('=')[1]) for x in splitsrc[2].split(',')})
        plot.Set(graphs[-1], **settings)
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), graphs[-1], True)
            DrawAxisHists(pads, axis, pads[0])
        graphs[-1].Draw('PZSAME')
        legend.AddEntry(graphs[-1], '', 'PL')


axis[0].GetYaxis().SetTitle('95%% CL limit on %s' % args.limit_on)
if args.y_title is not None:
    axis[0].GetYaxis().SetTitle(args.y_title)
axis[0].GetXaxis().SetTitle(args.x_title)
axis[0].GetXaxis().SetLabelOffset(axis[0].GetXaxis().GetLabelOffset()*2)

pads[0].SetLogy(1)

plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)
legend.Draw()

pads[0].cd()
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

# # CMS logo
plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)
# # plot.DrawTitle(pads[0], '0.1 fb^{-1} (13 TeV)', 3)

canv.Print('.png')
canv.Print('.pdf')



