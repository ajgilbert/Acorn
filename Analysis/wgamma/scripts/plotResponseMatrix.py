import ROOT
import argparse
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(True)
plot.ModTDRStyle(width=750, height=600, l=0.18, b=0.18, r=0.23)
# plot.SetDeepSeaPalette()
ROOT.gStyle.SetPalette(103)
ROOT.TColor.InvertPalette()
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(1.0)
ROOT.gStyle.SetPaintTextFormat('.2f')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', default='output_2016_eft_region_puppi_phi_f_binned.root')
parser.add_argument('--selection', default='eft_region')
parser.add_argument('--charge', default='p')
parser.add_argument('--channel', default='m')
parser.add_argument('--label', default='')
parser.add_argument('-o', '--output', default='response_matrix')
args = parser.parse_args()


f_in = ROOT.TFile(args.input)
hists = Node()
TDirToNode(f_in, '/', hists)

gen_ref = hists[args.channel]['XS']['2D']['WG_%s_%s_acc' % (args.charge, args.channel)]

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
                dirlabel = 'abs(reco_puppi_phi_f)' if args.selection == 'eft_region' else '0.5'
                h_dir = hists[args.channel]['%s_%s_%i' % (args.charge, args.channel, ix2)][dirlabel]
                if args.selection == 'eft_region':
                    h_title = 'WG_main_%s_%i_%i' % (args.charge, ix, iy)
                if args.selection == 'fid_region':
                    h_title = 'WG_main_%s_%i' % (args.charge, ix)
                if h_title in h_dir.d:
                    rec_events = h_dir[h_title].GetBinContent(iy2 + 1)
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
h_resp.GetYaxis().SetLabelSize(0.03)
h_resp.GetZaxis().SetLabelSize(0.03)
h_resp.GetZaxis().SetTitle('#varepsilon(reco./gen.)')
#
h_resp.SetContour(255)

if args.selection == 'eft_region':
    h_resp.GetXaxis().SetLabelOffset(h_resp.GetXaxis().GetLabelOffset() * 0.3)
    h_resp.GetYaxis().SetLabelOffset(h_resp.GetYaxis().GetLabelOffset() * 0.3)
    for ix in xrange(n_bins_x):
        h_resp.GetXaxis().SetBinLabel(ix * n_bins_y + 1, '[0, #pi/6]')
        h_resp.GetXaxis().SetBinLabel(ix * n_bins_y + 2, '[#pi/6, #pi/3]')
        h_resp.GetXaxis().SetBinLabel(ix * n_bins_y + 3, '[#pi/3, #pi/2]')
        h_resp.GetYaxis().SetBinLabel(ix * n_bins_y + 1, '[0, #pi/6]')
        h_resp.GetYaxis().SetBinLabel(ix * n_bins_y + 2, '[#pi/6, #pi/3]')
        h_resp.GetYaxis().SetBinLabel(ix * n_bins_y + 3, '[#pi/3, #pi/2]')
if args.selection == 'fid_region':
    h_resp.GetXaxis().SetLabelOffset(h_resp.GetXaxis().GetLabelOffset() * 0.5)
    h_resp.GetYaxis().SetLabelOffset(h_resp.GetYaxis().GetLabelOffset() * 0.5)
    h_resp.GetXaxis().SetLabelSize(0.04)
    h_resp.GetYaxis().SetLabelSize(0.04)
    bin_labels = ['[30, 50]', '[50,70]', '[70,100]', '[100,150]', '[150,200]', '[200,300]', '[300,500]', '[500,800]', '[800,1200]']
    for ib, label in enumerate(bin_labels):
        h_resp.GetXaxis().SetBinLabel(ib + 1, label)
        h_resp.GetYaxis().SetBinLabel(ib + 1, label)


h_resp.GetXaxis().LabelsOption('v')


canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()


ROOT.gStyle.SetTextFont(42)
h_resp.Draw('COLZ')


latex = ROOT.TLatex()

if args.selection == 'eft_region':
    x_offset = -2.5
    plot.Set(latex, TextSize=0.02, TextAlign=22)
    latex.DrawLatex(1, x_offset, '[150, 200]')
    latex.DrawLatex(4, x_offset, '[200, 300]')
    latex.DrawLatex(7, x_offset, '[300, 500]')
    latex.DrawLatex(10, x_offset, '[500, 800]')
    latex.DrawLatex(13, x_offset, '[800, 1200]')
    # latex.DrawLatex(16, x_offset, '[850, 1200]')

    plot.Set(latex, TextAngle=90)

    y_offset = -2.4
    latex.DrawLatex(y_offset, 1, '[150, 200]')
    latex.DrawLatex(y_offset, 4, '[200, 300]')
    latex.DrawLatex(y_offset, 7, '[300, 500]')
    latex.DrawLatex(y_offset, 10, '[500, 800]')
    latex.DrawLatex(y_offset, 13, '[800, 1200]')
    # latex.DrawLatex(y_offset, 16, '[850, 1200]')

plot.Set(latex, TextAngle=0, TextSize=0.03, TextAlign=21)
if args.selection == 'eft_region':
    latex.DrawLatex(12, -3.7, 'Gen. p_{T}^{#gamma} [GeV] x |#phi_{f}| bin')
if args.selection == 'fid_region':
    latex.DrawLatex(7.0, -2.4, 'Gen. p_{T}^{#gamma} [GeV] bin')
plot.Set(latex, TextAngle=90, TextSize=0.03, TextAlign=21)
if args.selection == 'eft_region':
    latex.DrawLatex(-3.2, 12, 'Reco. p_{T}^{#gamma} [GeV] x |#phi_{f}| bin')
if args.selection == 'fid_region':
    latex.DrawLatex(-2.0, 7.0, 'Reco. p_{T}^{#gamma} [GeV] bin')

plot.FixOverlay()

plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.12, 0.035, 1.2)
plot.DrawTitle(pads[0], args.label, 3)


canv.Print('.pdf')
canv.Print('.png')

