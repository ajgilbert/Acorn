import ROOT
import argparse
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(True)
plot.ModTDRStyle(width=750, height=600, l = 0.18, b = 0.18, r = 0.23)
# plot.SetDeepSeaPalette()
ROOT.gStyle.SetPalette(103)
ROOT.TColor.InvertPalette()
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(1.0)
ROOT.gStyle.SetPaintTextFormat('.2f')

# parser = argparse.ArgumentParser()
# parser.add_argument('-i', help="input.root:tree")
# parser.add_argument('-d', default='', help='Draw expression')
# parser.add_argument('-b', default='', help='Binning')
# parser.add_argument('-s', default='', help='Selection/weights')
# parser.add_argument('-o', default='', help='output.root:hist')
# args = parser.parse_args()


f_in = ROOT.TFile('output_2016_eft_region_puppi_phi_f_binned.root')
hists = Node()
TDirToNode(f_in, '/', hists)

gen_ref = hists['m']['XS']['2D']['WG_p_m_acc']

gen_ref.Print('range')

n_bins_x = gen_ref.GetNbinsX()
n_bins_y = gen_ref.GetNbinsY()
n_bins = n_bins_x * n_bins_y

h_resp  = ROOT.TH2F('response_matrix', '', n_bins, -0.5, float(n_bins) - 0.5, n_bins, -0.5, float(n_bins) - 0.5)

for iy in xrange(n_bins_y):
    for ix in xrange(n_bins_x):
        gen_events = gen_ref.GetBinContent(ix + 1, iy + 1)
        for iy2 in xrange(n_bins_y):
            for ix2 in xrange(n_bins_x):
                # Does the rec hist even exist?
                h_dir = hists['m']['p_m_%i' % ix2]['abs(reco_puppi_phi_f)']
                h_title = 'WG_main_p_%i_%i' % (ix, iy)
                if h_title in h_dir.d:
                    rec_events = h_dir['WG_main_p_%i_%i' % (ix, iy)].GetBinContent(iy2 + 1)
                    if rec_events < 0.:
                        rec_events = 0.
                else:
                    rec_events = 0
                h_resp.SetBinContent(1 + iy + ix * n_bins_y, 1 + iy2 + ix2 * n_bins_y, rec_events / gen_events)

h_resp.Print('range')


h_resp.GetXaxis().SetLabelFont(42)
h_resp.GetXaxis().SetTickLength(0)
h_resp.GetYaxis().SetTickLength(0)
h_resp.GetYaxis().SetLabelFont(42)
h_resp.GetXaxis().SetLabelSize(0.03)
h_resp.GetXaxis().SetLabelOffset(h_resp.GetXaxis().GetLabelOffset() * 0.3)
h_resp.GetYaxis().SetLabelOffset(h_resp.GetYaxis().GetLabelOffset() * 0.3)
h_resp.GetYaxis().SetLabelSize(0.03)
h_resp.GetZaxis().SetLabelSize(0.03)
# h_resp.GetZaxis().SetTitle('N_{reco}^{j} / N_{true}^{i}')
h_resp.GetZaxis().SetTitle('#varepsilon(reco./gen.)')
#
h_resp.SetContour(255)

for ix in xrange(n_bins_x):
    h_resp.GetXaxis().SetBinLabel(ix * n_bins_y + 1, '[0, #pi/6]')
    h_resp.GetXaxis().SetBinLabel(ix * n_bins_y + 2, '[#pi/6], #pi/3]')
    h_resp.GetXaxis().SetBinLabel(ix * n_bins_y + 3, '[#pi/3], #pi/2]')
    h_resp.GetYaxis().SetBinLabel(ix * n_bins_y + 1, '[0, #pi/6]')
    h_resp.GetYaxis().SetBinLabel(ix * n_bins_y + 2, '[#pi/6], #pi/3]')
    h_resp.GetYaxis().SetBinLabel(ix * n_bins_y + 3, '[#pi/3], #pi/2]')

h_resp.GetXaxis().LabelsOption('v')


canv = ROOT.TCanvas('response_matrix', 'response_matrix')
pads = plot.OnePad()


ROOT.gStyle.SetTextFont(42)
# hist.SetMarkerSize(args.marker_size)
h_resp.Draw('COLZ')


latex = ROOT.TLatex()
x_offset = -3.0
plot.Set(latex, TextSize=0.02, TextAlign=22)
latex.DrawLatex(1, x_offset, '[150, 210]')
latex.DrawLatex(4, x_offset, '[210, 300]')
latex.DrawLatex(7, x_offset, '[300, 420]')
latex.DrawLatex(10, x_offset, '[420, 600]')
latex.DrawLatex(13, x_offset, '[600, 850]')
latex.DrawLatex(16, x_offset, '[850, 1200]')
# latex.DrawLatex(8.5, -2.0, 'VBF')
# latex.DrawLatex(13, -2.0, 'WH')
# latex.DrawLatex(17, -2.0, 'ZH')
# latex.DrawLatex(21.5, -2.0, 'ttH')
plot.Set(latex, TextAngle=90)

y_offset = -2.6
latex.DrawLatex(y_offset, 1, '[150, 210]')
latex.DrawLatex(y_offset, 4, '[210, 300]')
latex.DrawLatex(y_offset, 7, '[300, 420]')
latex.DrawLatex(y_offset, 10, '[420, 600]')
latex.DrawLatex(y_offset, 13, '[600, 850]')
latex.DrawLatex(y_offset, 16, '[850, 1200]')
# latex.DrawLatex(-2.5, 2.5, 'ttH')
# latex.DrawLatex(-2.5, 7, 'ZH')
# latex.DrawLatex(-2.5, 11, 'WH')
# latex.DrawLatex(-2.5, 15.5, 'VBF')
# latex.DrawLatex(-2.5, 21, 'ggH')
plot.Set(latex, TextAngle=0, TextSize=0.03, TextAlign=21)
latex.DrawLatex(15, -4.2, 'Gen. p_{T}^{#gamma} [GeV] x |#phi_{f}| bin')
plot.Set(latex, TextAngle=90, TextSize=0.03, TextAlign=21)
latex.DrawLatex(-3.7, 15, 'Reco. p_{T}^{#gamma} [GeV] x |#phi_{f}| bin')

plot.FixOverlay()

plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.12, 0.035, 1.2)
# if args.label is not None: plot.DrawTitle(pads[0], args.label, 3)


canv.Print('.pdf')
canv.Print('.png')
# fin = ROOT.TFile(args.i.split(':')[0])
# tree = fin.Get(args.i.split(':')[1])

# tree.Draw('%s>>h%s' % (args.d, args.b), args.s)

# h = ROOT.gDirectory.Get('h')
# h.SetName(args.o.split(':')[1])

# fout = ROOT.TFile(args.o.split(':')[0], 'RECREATE')
# h.Write()
# fout.Close()

# fin.Close()
