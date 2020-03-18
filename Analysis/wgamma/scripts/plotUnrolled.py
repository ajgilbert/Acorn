import ROOT
import argparse
import json
import itertools
from copy import deepcopy
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
from Acorn.Analysis.plottemplates import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)

parser = argparse.ArgumentParser()
parser.add_argument('--selection', default='eft_region')
parser.add_argument('--channel', default='l')
parser.add_argument('--charge', default='x')
parser.add_argument('--year', default='total')
parser.add_argument('--output', '-o', default='postfit_plot')
# parser.add_argument('--scheme', default='phi_f_binned')
# parser.add_argument('--years', default='2016,2017,2018')
# parser.add_argument('--ratio', action='store_true')
# parser.add_argument('--charge', default='p')
args = parser.parse_args()

lumi_year = {
    '2016': '35.9',
    '2017': '41.5',
    '2018': '59.3',
    'total': '136.9'
}

if args.selection == 'eft_region':
    plot.ModTDRStyle(width=1000)
if args.selection == 'fid_region':
    plot.ModTDRStyle(width=600)

ROOT.gStyle.SetEndErrorSize(0)

def DoScaleEnvelope(node, nominal):
    h = node[nominal].Clone()
    h_alts = []
    for i in range(6):
        h_alts.append(node[nominal + '__sc_%i' % i])
    for ix in xrange(1, h.GetNbinsX() + 1):
        for iy in xrange(1, h.GetNbinsY() + 1):
            max_dev = max([abs(x.GetBinContent(ix, iy) - h.GetBinContent(ix, iy)) for x in h_alts])
            h.SetBinError(ix, iy, max_dev)
    return h

def DivideGraphByHist(gr, hist):
    res = gr.Clone()
    for i in xrange(gr.GetN()):
        res.GetY()[i] = res.GetY()[i] / hist.GetBinContent(i + 1)
    if type(res) is ROOT.TGraphAsymmErrors:
        for i in xrange(gr.GetN()):
            res.GetEYhigh()[i] = res.GetEYhigh()[i]/hist.GetBinContent(i + 1)
            res.GetEYlow()[i] = res.GetEYlow()[i]/hist.GetBinContent(i + 1)
    return res


def RebinHists(h_names, h_fits, pt_bins, phi_bins, chg, chn, year, fit):
    h_dicts = []
    n_bins_phi = len(phi_bins) - 1
    n_bins_pt = len(pt_bins) - 1
    for i in xrange(n_bins_phi):
        h_dict = {}
        for h_name in h_names:
            h_dict[h_name] = ROOT.TH1F(h_name, h_name, len(pt_bins) - 1, array('d', pt_bins))
            for ib in xrange(n_bins_pt):
                src_dir = '%s_%s_%i_%s_%s' % (chg, chn, ib, year, fit)
                # It might be that the proc yield was zero, so never went into the cards
                if h_name in h_fits[src_dir].d:
                    h_dict[h_name].SetBinContent(ib + 1, h_fits[src_dir][h_name].GetBinContent(i + 1))
                    h_dict[h_name].SetBinError(ib + 1, h_fits[src_dir][h_name].GetBinError(i + 1))
                else:
                    for dname in h_fits[src_dir].d:
                        if dname.startswith('WG_') and h_name in dname:
                            h_dict[h_name].SetBinContent(ib + 1, h_dict[h_name].GetBinContent(ib + 1) + h_fits[src_dir][dname].GetBinContent(i + 1))
                            curr_err_2 = math.pow(h_dict[h_name].GetBinContent(ib + 1), 2)
                            h_dict[h_name].SetBinError(ib + 1, math.sqrt(curr_err_2 + math.pow(h_fits[src_dir][dname].GetBinError(i + 1), 2)))

        h_dicts.append(h_dict)
    return h_dicts


chg = args.charge
chn = args.channel

if args.selection == 'eft_region':
    f_fits = ROOT.TFile('shapes_combined_puppi_phi_f_binned.root')
    h_fits = Node()
    TDirToNode(f_fits, '/', h_fits)

    f_fits_c3w_p1 = ROOT.TFile('shapes_combined_puppi_phi_f_binned_c3w_p1.root')
    h_fits_c3w_p1 = Node()
    TDirToNode(f_fits_c3w_p1, '/', h_fits_c3w_p1)

    f_fits_c3w_m1 = ROOT.TFile('shapes_combined_puppi_phi_f_binned_c3w_m1.root')
    h_fits_c3w_m1 = Node()
    TDirToNode(f_fits_c3w_m1, '/', h_fits_c3w_m1)

    year = args.year
    fit = 'prefit'

    pt_bins = BinEdgesFromStr('[150,200,300,500,800,1200]')
    phi_bins = BinEdgesFromStr('(3,0.,math.pi/2.)')

    h_names = ["TotalProcs", "data_obs", "VV_R", "VV_E", "DY_XZG_R", "ZG_IZG_R", "DY_E", "TTG_ITTG_R", "TT_XTTG_R", "TT_E" , "data_fakes_highpt" , "GG_R", "GG_E" , "WG_ooa_p", "WG_ooa_n" , "WG_met1_p", "WG_met1_n" , "WG_main_p", "WG_main_n"]

if args.selection == 'fid_region':
    f_fits = ROOT.TFile('shapes_combined_pt_diff_fid_pt_binned.root')
    h_fits = Node()
    TDirToNode(f_fits, '/', h_fits)
    year = args.year
    fit = 'prefit'

    pt_bins = BinEdgesFromStr('[30,50,70,100,150,200,300,500,800,1200]')
    phi_bins = BinEdgesFromStr('(1,0.,1.)')

    h_names = ["TotalProcs", "data_obs", "VV_R", "VV_E", "DY_XZG_R", "ZG_IZG_R", "DY_E", "TTG_ITTG_R", "TT_XTTG_R", "TT_E" , "data_fakes_sub" , "data_fakes_lep_sub", "GG_R", "GG_E" , "WG_ooa_x", "WG_met1_x" , "WG_main_x"]



chg_latex = {
    'n': '-',
    'p': '+',
    'x': '#pm'
}
chn_latex = {
    'm': '#mu',
    'e': 'e',
    'l': 'l'
}
n_bins_pt = len(pt_bins) - 1
n_bins_phi = len(phi_bins) - 1
print pt_bins, phi_bins
# h_names = h_fits['%s_%s_0_%s_%s' % (chg, chn, year, fit)].d.keys()

h_dicts = RebinHists(h_names, h_fits, pt_bins, phi_bins, chg, chn, year, fit)

if args.selection == 'eft_region':
    h_dicts_c3w_p1 = RebinHists(h_names, h_fits_c3w_p1, pt_bins, phi_bins, chg, chn, year, fit)
    h_dicts_c3w_m1 = RebinHists(h_names, h_fits_c3w_m1, pt_bins, phi_bins, chg, chn, year, fit)

    for i in xrange(len(h_dicts)):
        h_dicts[i]['TotalProcs_c3w_p1'] = h_dicts_c3w_p1[i]['TotalProcs']
        h_dicts[i]['TotalProcs_c3w_m1'] = h_dicts_c3w_m1[i]['TotalProcs']

width = 1. / n_bins_phi

canv = ROOT.TCanvas(args.output, args.output)

ratio = None
# if args.ratio:
#     ratio = 0.4
pads, ratio_pads, purity_pads = SetupPads([width] * (n_bins_phi - 1), [0, 0], [0, 0], ratio=0.20, purity=0.20)


with open('input/wgamma_plot_layouts.json') as jsonfile:
    layouts = json.load(jsonfile)

plotcfg = DEFAULT_CFG
plotcfg.update({
    'type': 'datamc',
    'ratio': True,
    'purity': True,
    'global_hist_opts': {
        'draw_opts': 'HIST',
        'legend_opts': 'F',
        'marker_size': 0.6,
        'line_width': 1
    },
    'legend_pos': [0.80, 0.60, 0.95, 0.88],
    'main_logo': '',
    'sub_logo': '',
    'top_title_right': '',
    'ratio_y_range': [0.31, 1.79]
})
if args.selection == 'fid_region':
    plotcfg['legend_pos'] = [0.65, 0.62, 0.95, 0.93]
    plotcfg['ratio_y_range'] = [0.61, 1.39]

print plotcfg

out_objects = []

latex = ROOT.TLatex()
latex.SetTextFont(62)
latex.SetTextSize(0.03)
latex.SetTextAlign(22)

for i in xrange(n_bins_phi):
    pads[i].cd()

    if args.selection == 'eft_region':
        thiscfg = UpdateDict(plotcfg, {
                              'pads': [pads[i], purity_pads[i], ratio_pads[i]],
                              'overlays': [
                                {"name": "c3wUp", "entries": ["TotalProcs_c3w_p1"], 'hist_postfix': '', 'legend': 'C_{3W} = 0.2 TeV^{-2}', 'color': 4},
                                {"name": "c3wDown", "entries": ["TotalProcs_c3w_m1"], 'hist_postfix': '', 'legend': 'C_{3W} = -0.2  TeV^{-2}', 'color': 2}
                              ]
                              })
        layout_name = 'data_fakes_EFT_simple'
    else:
        thiscfg = UpdateDict(plotcfg, {
                              'pads': [pads[i], purity_pads[i], ratio_pads[i]],
                              })
        layout_name = 'data_fakes_diff_simple'

    if i == 0:
        thiscfg['main_logo'] = 'CMS'
        thiscfg['sub_logo'] = 'Internal'
        thiscfg['top_title_left'] = 'W^{%s}(%s^{%s}#nu)#gamma' % (chg_latex[chg], chn_latex[chn], chg_latex[chg])
    if i == n_bins_phi - 1:
        thiscfg['top_title_right'] = '%s fb^{-1} (13 TeV)' % lumi_year[year]
    h_dicts[i]['TotalProcs'].Print("range")
    res = MakeMultiHistPlot('test',
                      outdir='.',
                      hists=h_dicts[i],
                      cfg=thiscfg,
                      layout=layouts[layout_name],
    )

    out_objects.append(res)
    # h_dicts_c3w_p1[i]['TotalProcs'].Scale(1.0, 'width')
    # plot.Set(h_dicts_c3w_p1[i]['TotalProcs'], LineColor=4, LineWidth=3)
    # h_dicts_c3w_p1[i]['TotalProcs'].Draw('HISTSAME')

    pads[i].cd()
    pads[i].SetLogy(True)

    h = plot.GetAxisHist(pads[i])
    hr = plot.GetAxisHist(ratio_pads[i])
    if args.selection == 'eft_region':
        h.SetMinimum(1E-4)
        h.SetMaximum(1E+2)
    if args.selection == 'fid_region':
        h.SetMinimum(1E-2)
        h.SetMaximum(1E+5)

    # h.Draw()
    h.GetXaxis().SetNdivisions(510)
    hr.GetXaxis().SetNdivisions(510)
    hr.GetXaxis().ChangeLabel(-1, -1., -1., -1, -1, -1, ' ')
    h.GetYaxis().SetTickLength(h.GetYaxis().GetTickLength() * 0.5)
    hr.GetYaxis().SetTickLength(hr.GetYaxis().GetTickLength() * 0.5)
    if i == 0:
        h.GetYaxis().SetTitle('Events / GeV')
    if i > 0:
        h.GetYaxis().SetLabelSize(0)
        hr.GetYaxis().SetLabelSize(0)
        h.GetYaxis().SetTitle('')
        hr.GetYaxis().SetTitle('')
    if i == n_bins_phi - 1:
        hr.GetXaxis().SetTitle('Photon p_{T} (GeV)')


    pad_width = 1. - pads[i].GetLeftMargin() - pads[i].GetRightMargin()
    if args.selection == 'eft_region':
        latex.DrawLatexNDC(pads[i].GetLeftMargin() + pad_width * 0.75, 0.91, '%.2f #leq |#phi_{f}| < %.2f' % (phi_bins[i], phi_bins[i + 1]))


canv.Print('.png')
canv.Print('.pdf')
