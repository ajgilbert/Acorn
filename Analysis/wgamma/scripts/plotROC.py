#!/usr/bin/env python
import ROOT
import CombineHarvester.CombineTools.plotting as plot
import argparse
from array import array

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)


def ColorLabels(text):
    return text


parser = argparse.ArgumentParser()
# parser.add_argument('--output', '-o', help='output name')
# parser.add_argument('--input', '-i', help='Main input file')
# parser.add_argument('--label', '-l', help='Main input file')
# parser.add_argument('--label-size', type=float, default=0.04)
# parser.add_argument('--marker-size', type=float, default=1.5)
# parser.add_argument('--print-twiki', action='store_true')
# parser.add_argument('--titles', default='x,y,z')
# parser.add_argument('--draw-opt', default='COLZTEXT')
# parser.add_argument('--palette', default='default', choices=['default', 'i103'])
# parser.add_argument(
#     '--subline', default='Internal', help='text to add next to cms logo')

args = parser.parse_args()
plot.ModTDRStyle()

fin = ROOT.TFile('output_2018_photon_fakes_200511.root')
# barrel_m_iso_t_sig_t_pt_30_40

def MakeCurve(file, barrel, cat, maxiso=1000.):
    h_WG = fin.Get('/m/%s_m_iso_t_sig_t_pt_%s/p0_worstiso/WG_R' % (barrel, cat))
    h_WJ = fin.Get('/m/%s_m_iso_t_sig_t_pt_%s/p0_worstiso/W_J' % (barrel, cat))
    h_WG_e = fin.Get('/e/%s_e_iso_t_sig_t_pt_%s/p0_worstiso/WG_R' % (barrel, cat))
    h_WJ_e = fin.Get('/e/%s_e_iso_t_sig_t_pt_%s/p0_worstiso/W_J' % (barrel, cat))
    h_WG.Add(h_WG_e)
    h_WJ.Add(h_WJ_e)
    i_WG = h_WG.Integral(0, h_WG.GetNbinsX() + 1)
    i_WJ = h_WJ.Integral(0, h_WJ.GetNbinsX() + 1)

    xvals = []
    yvals = []

    for i in xrange(1, h_WG.GetNbinsX() + 1):
        edge = h_WG.GetXaxis().GetBinUpEdge(i)
        if edge > maxiso:
            continue
        eff_WG = h_WG.Integral(0, i) / i_WG
        rej_WJ = h_WJ.Integral(i + 1, h_WJ.GetNbinsX() + 1) / i_WJ

        xvals.append(eff_WG)
        yvals.append(rej_WJ)
        print '>> %f, %f, %f' % (edge, eff_WG, rej_WJ)

    gr = ROOT.TGraph(len(xvals), array('d', xvals), array('d', yvals))
    return gr


# plot.ModTDRStyle(width=750, height=600, l = 0.17, b = 0.15, r = 0.23, t=0.08)

# plot.SetCorrMatrixPalette()

# ROOT.gStyle.SetNdivisions(510, "XYZ")
# ROOT.gStyle.SetMarkerSize(1.0)
# ROOT.gStyle.SetPaintTextFormat('.2f')

# x_title, y_title, z_title = args.titles.split(',')

# inputs = args.input.split(':')
# fin = ROOT.TFile(inputs[0])
# hist = fin.Get(inputs[1])

# hist.GetXaxis().SetTitle(x_title)
# hist.GetYaxis().SetTitle(y_title)
# hist.GetZaxis().SetTitle(z_title)
# hist.GetXaxis().SetLabelFont(42)
# hist.GetYaxis().SetLabelFont(42)
# hist.GetXaxis().SetTickLength(0)
# hist.GetYaxis().SetTickLength(0)
# hist.GetXaxis().SetLabelSize(args.label_size)
# hist.GetYaxis().SetLabelSize(args.label_size)
# hist.GetZaxis().SetLabelSize(0.045)
# # hist.GetZaxis().SetTitle('#rho')
# hist.SetContour(255)

canv = ROOT.TCanvas('worstiso_roc', 'worstiso_roc')
pads = plot.OnePad()


targets = [
    ('30_40', 1, '30-40 GeV'),
    ('40_60', 2, '40-60 GeV'),
    ('60_100', 4, '60-100 GeV'),
    ('100_200', 8, '100-200 GeV'),
]

legend = ROOT.TLegend(0.6, 0.86 - 0.04 * 3, 0.90, 0.91, '', 'NBNDC')

grs = []

for t in targets:
    grs.append(MakeCurve(fin, 'barrel', t[0], 12.0))
    plot.Set(grs[-1], LineColor=t[1], MarkerColor=t[1])
    legend.AddEntry(grs[-1], t[2], 'LP')

grs[0].Draw('APL')
for gr in grs:
    gr.Draw('LPSAME')

legend.Draw()
# ROOT.gStyle.SetTextFont(42)
# hist.SetMarkerSize(args.marker_size)
# hist.Draw(args.draw_opt)

# plot.FixOverlay()

plot.DrawCMSLogo(pads[0], 'CMS', 'Internal', 0, 0.17, 0.035, 1.2)

# if args.label is not None:
#     plot.DrawTitle(pads[0], args.label, 3)


canv.Print('.pdf')
canv.Print('.png')
