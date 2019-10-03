import ROOT
import argparse
import json
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
parser.add_argument('--scheme', default='phi_f_binned')
parser.add_argument('--years', default='2016,2017,2018')
parser.add_argument('--ratio', action='store_true')
parser.add_argument('--charge', default='p')
args = parser.parse_args()

settings = {
    'eft_region': {
        'canvas_width': 1000,
        'y_axis_title': '#Delta^{2}#sigma/#Deltap_{T}^{#gamma}#Delta|#phi_{f}| (fb/GeV)',
        'show_matrix': False,
        'y_axis_min': 5E-5,
        'y_axis_max': 8E-1
    },
    'fid_region': {
        'canvas_width': 600,
        'y_axis_title': '#Delta#sigma/#Deltap_{T}^{#gamma} (fb/GeV)',
        'show_matrix': True,
        'y_axis_min': 1E-1,
        'y_axis_max': 1E+4
    }
}


plot.ModTDRStyle(width=settings[args.selection]['canvas_width'])
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

def ZeroErrors(hist):
    for i in xrange(1, hist.GetNbinsX() + 1):
        hist.SetBinError(i, 1E-20)
    return hist

years = args.years.split(',')
chg = args.charge

ref_hists = {}

# These factors undo the kNNLO that is applied for the nominal prediction
corr_factors = {
    # '2016': 178.9 / 214.68,
    # '2017': 191.4 / 229.92,
    # '2018': 191.4 / 229.92
    '2016': 1.0,
    '2017': 1.0,
    '2018': 1.0
}

doScaleUncert = True

for yr in years:
    f_ref = ROOT.TFile('output_%s_%s_%s.root' % (yr, args.selection, args.scheme))
    h_ref_m = Node()
    TDirToNode(f_ref, 'm/XS/2D', h_ref_m)
    h_ref_e = Node()
    TDirToNode(f_ref, 'e/XS/2D', h_ref_e)

    if doScaleUncert:
        h_xs_m = DoScaleEnvelope(h_ref_m, 'XS_WG_%s_m_acc' % chg)
        h_xs_e = DoScaleEnvelope(h_ref_e, 'XS_WG_%s_e_acc' % chg)
    else:
        h_xs_m = h_ref_m['XS_WG_%s_m_acc' % chg]
        h_xs_e = h_ref_e['XS_WG_%s_e_acc' % chg]

    h_xs = h_xs_m.Clone()
    h_xs.Add(h_xs_e)

    ## Assume these uncertainties are fully correlated - replace the
    ## automatic error with the linear sum
    for ix in xrange(1, h_xs.GetNbinsX() + 1):
        for iy in xrange(1, h_xs.GetNbinsY() + 1):
            h_xs.SetBinError(ix, iy, h_xs_m.GetBinError(ix, iy) + h_xs_e.GetBinError(ix, iy))


    # Scale from pb to fb, then a factor 3/2 since we only added e + mu, and we want all l
    h_xs.Scale(1000. * (3. / 2.))
    h_xs.Scale(corr_factors[yr])
    h_xs.Print("range")
    h_xs.Scale(1, 'width')
    h_xs.Print("range")

    ref_hists[yr] = h_xs

n_bins_pt = ref_hists[years[0]].GetNbinsX()
n_bins_phi = ref_hists[years[0]].GetNbinsY()
width = 1. / n_bins_phi

canv = ROOT.TCanvas('diff_xsec', 'diff_xsec')

ratio = None
if args.ratio:
    ratio = 0.2
pads, ratio_pads = SetupPads([width] * (n_bins_phi - 1), [0, 0], [0, 0], ratio=ratio)

xsec_results = {}
with open('%s.json' % args.scheme) as jsonfile:
    xsec_results = json.load(jsonfile)['xsec2D']

ref_hists_1D = []
obs_graphs = []

for i in range(n_bins_phi):
    ref_hists_1D.append(ref_hists[yr].ProjectionX('proj_%i' % i, i + 1, i + 1, 'e'))
    h_1D = ref_hists_1D[-1]
    x_vals = []
    y_vals = []
    ex_hi = []
    ex_lo = []
    ey_hi = []
    ey_lo = []

    for j in range(n_bins_pt):
        x_vals.append(h_1D.GetXaxis().GetBinCenter(j + 1))
        # ex_lo.append(h_1D.GetXaxis().GetBinWidth(j + 1) / 2.)
        # ex_hi.append(h_1D.GetXaxis().GetBinWidth(j + 1) / 2.)
        ex_lo.append(0.)
        ex_hi.append(0.)
        if args.selection == 'fid_region':
            poi_res = xsec_results['r_%s_%i' % (chg, j)]
        else:
            poi_res = xsec_results['r_%s_%i_%i' % (chg, j, i)]
        # poi_res = xsec_results['r_%s_%i' % (chg, j)]
        y_vals.append(h_1D.GetBinContent(j + 1) * poi_res['Val'])
        ey_hi.append(h_1D.GetBinContent(j + 1) * poi_res['ErrorHi'])
        ey_lo.append(h_1D.GetBinContent(j + 1) * poi_res['ErrorLo'] * -1.)
    obs_graphs.append(ROOT.TGraphAsymmErrors(len(x_vals), array('d', x_vals), array('d', y_vals), array('d', ex_lo), array('d', ex_hi), array('d', ey_lo), array('d', ey_hi)))

h_axes = [h.Clone() for h in ref_hists_1D]
r_h_axes = [h.Clone() for h in ref_hists_1D]

legend = ROOT.TLegend(0.68, 0.72, 0.94, 0.88, '', 'NBNDC')

h_obs = ROOT.TH1F('h_obs', '', 1, 0, 1)
plot.Set(h_obs, LineWidth=2)

latex = ROOT.TLatex()
latex.SetTextFont(62)
latex.SetTextSize(0.03)
latex.SetTextAlign(22)
# latex.SetTextColor(14)
h_matrix = None
h_matrix_fill = None

h_store = {}

for i, h in enumerate(h_axes):
    print i
    hr = r_h_axes[i]
    h.Reset()
    hr.Reset()

    pads[i].cd()
    pads[i].SetLogy(True)
    h.Draw()

    ratio_pads[i].cd()
    hr.Draw()

    pads[i].cd()
    plot.SetupTwoPadSplitAsRatio(
        [pads[i], ratio_pads[i]], h, hr, '(Obs-Exp)/Exp', True, -0.28, +0.28)

    h.SetMinimum(settings[args.selection]['y_axis_min'])
    h.SetMaximum(settings[args.selection]['y_axis_max'])

    h.GetXaxis().SetNdivisions(510)
    hr.GetXaxis().SetNdivisions(510)
    hr.GetXaxis().ChangeLabel(-1, -1., -1., -1, -1, -1, ' ')
    h.GetYaxis().SetTickLength(h.GetYaxis().GetTickLength() * 0.5)
    hr.GetYaxis().SetTickLength(hr.GetYaxis().GetTickLength() * 0.5)
    if i == 0:
        h.GetYaxis().SetTitle(settings[args.selection]['y_axis_title'])
    if i > 0:
        h.GetYaxis().SetLabelSize(0)
        hr.GetYaxis().SetLabelSize(0)
        h.GetYaxis().SetTitle('')
        hr.GetYaxis().SetTitle('')
    if i == n_bins_phi - 1:
        hr.GetXaxis().SetTitle('Photon p_{T} (GeV)')

    plot.Set(ref_hists_1D[i], LineWidth=1, LineColor=2, MarkerSize=0, FillColor=plot.CreateTransparentColor(2, 0.2))
    h_store['ref_hists_1D_%i_line' % i] = ZeroErrors(ref_hists_1D[i].Clone())
    plot.Set(obs_graphs[i], LineWidth=2, MarkerSize=0.6)
    ref_hists_1D[i].Draw('E2SAME')
    h_store['ref_hists_1D_%i_line' % i].Draw('LSAME')


    if i == 0:
        legend.AddEntry(h_obs, 'Observed', 'PE')
        legend.AddEntry(ref_hists_1D[i], 'aMC@NLO', 'LF')


    ##### THIS PART ONLY FOR THE diff XSEC
    if settings[args.selection]['show_matrix']:
        pt_binning = []
        for ib in xrange(1, ref_hists_1D[i].GetNbinsX() + 1):
            pt_binning.append(ref_hists_1D[i].GetXaxis().GetBinLowEdge(ib))
            if ib == ref_hists_1D[i].GetNbinsX():
                pt_binning.append(ref_hists_1D[i].GetXaxis().GetBinUpEdge(ib))
        matrix_p = "/afs/cern.ch/work/a/agilbert/matrix/MATRIX_v1.0.3/run/ppexnea03_MATRIX/result/run_NNLO_iso_met/gnuplot/histograms/pT_gamma__NNLO_QCD.hist"
        matrix_m = "/afs/cern.ch/work/a/agilbert/matrix/MATRIX_v1.0.3/run/ppenexa03_MATRIX/result/run_NNLO_iso_met/gnuplot/histograms/pT_gamma__NNLO_QCD.hist"
        h_matrix_p = ReadTxtHist(matrix_p)
        h_matrix_m = ReadTxtHist(matrix_m)
        h_matrix_p_sc_lo = ReadTxtHist(matrix_p, column=3)
        h_matrix_p_sc_hi = ReadTxtHist(matrix_p, column=5)
        h_matrix_m_sc_lo = ReadTxtHist(matrix_m, column=3)
        h_matrix_m_sc_hi = ReadTxtHist(matrix_m, column=5)
        h_matrix = h_matrix_p.Clone()
        h_matrix_sc_lo = h_matrix_p_sc_lo.Clone()
        h_matrix_sc_hi = h_matrix_p_sc_hi.Clone()
        h_matrix.Add(h_matrix_m)
        h_matrix_sc_lo.Add(h_matrix_m_sc_lo)
        h_matrix_sc_hi.Add(h_matrix_m_sc_hi)

        h_matrix.Scale(3)
        h_matrix = VariableRebin(h_matrix, pt_binning)
        h_matrix_sc_lo.Scale(3)
        h_matrix_sc_lo = VariableRebin(h_matrix_sc_lo, pt_binning)
        h_matrix_sc_hi.Scale(3)
        h_matrix_sc_hi = VariableRebin(h_matrix_sc_hi, pt_binning)

        for ib in xrange(1, h_matrix.GetNbinsX() + 1):
            h_matrix.SetBinError(ib, max(abs(h_matrix.GetBinContent(ib) - h_matrix_sc_lo.GetBinContent(ib)), abs(h_matrix.GetBinContent(ib) - h_matrix_sc_hi.GetBinContent(ib))))
        h_matrix.Scale(1., 'width')
        print h_matrix.Integral(), ref_hists_1D[i].Integral()
        plot.Set(h_matrix, LineColor=4, LineWidth=1, MarkerSize=0, FillColorAlpha=(4, 0.3))
        h_store['matrix_%i_line' % i] = ZeroErrors(h_matrix.Clone())
        h_matrix.Draw('E2SAME')
        h_store['matrix_%i_line' % i].Draw('LSAME')
        legend.AddEntry(h_matrix, 'MATRIX (NNLO QCD)', 'LF')


    obs_graphs[i].Draw('SAMEP')


    pad_width = 1. - pads[i].GetLeftMargin() - pads[i].GetRightMargin()
    if args.selection != 'fid_region':
        latex.DrawLatexNDC(pads[i].GetLeftMargin() + pad_width * 0.75, 0.91, '%.2f #leq |#phi_{f}| < %.2f' % (ref_hists[yr].GetYaxis().GetBinLowEdge(i + 1), ref_hists[yr].GetYaxis().GetBinUpEdge(i + 1)))

r_ref_hists_1D = []
r_obs_graphs = []
r_matrix = None
if args.ratio:
    for i, h in enumerate(r_h_axes):
        ratio_pads[i].cd()
        ratio_pads[i].SetGrid(0, 1)


        r_ref_hists_1D.append(plot.MakeRatioHist(ref_hists_1D[i], ref_hists_1D[i], True, False))
        r_obs_graphs.append(DivideGraphByHist(obs_graphs[i], ref_hists_1D[i]))
        for ib in xrange(1, r_ref_hists_1D[i].GetNbinsX() + 1):
            r_ref_hists_1D[i].SetBinContent(ib, r_ref_hists_1D[i].GetBinContent(ib) - 1.0)
        for ib in xrange(r_obs_graphs[i].GetN()):
            r_obs_graphs[i].GetY()[ib] = r_obs_graphs[i].GetY()[ib] - 1.0

        r_ref_hists_1D[i].Draw('E2SAME')
        h_store['r_ref_hists_1D_%i_line' % i] = ZeroErrors(r_ref_hists_1D[i].Clone())
        h_store['r_ref_hists_1D_%i_line' % i].Draw('LSAME')
        if settings[args.selection]['show_matrix']:
            h_store['r_matrix'] = plot.MakeRatioHist(h_matrix, ref_hists_1D[i], True, False)
            for ib in xrange(1, h_store['r_matrix'].GetNbinsX() + 1):
                h_store['r_matrix'].SetBinContent(ib, h_store['r_matrix'].GetBinContent(ib) - 1.0)
            h_store['r_matrix_line'] = ZeroErrors(h_store['r_matrix'].Clone())
            h_store['r_matrix'].Draw('E2SAME')
            h_store['r_matrix_line'].Draw('LSAME')

        r_obs_graphs[i].Draw('SAMEP')
        # r_ref_hists_1D[i].Draw('LSAME')

latex.SetTextSize(0.05)
chg_labels = {
    'p': '+',
    'n': '-',
    'x': '#pm'
}
# latex.DrawLatexNDC(0.17, 0.17, 'W^{%s}(#rightarrowl^{%s}#nu)#gamma' % (chg_labels[args.charge], chg_labels[args.charge]))

legend.Draw()


pads[0].cd()
# plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.24, 0.035, 1.2, cmsTextSize=0.9)
plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 11, 0.045, 0.05, 1.0, '', 1.0)

plot.DrawTitle(pads[-1], '136.9 fb^{-1} (13 TeV)', 3)
plot.DrawTitle(pads[0], 'W^{%s}(l^{%s}#nu)#gamma' % (chg_labels[args.charge], chg_labels[args.charge]), 1)

canv.Print('.png')
canv.Print('.pdf')
