import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from Acorn.Analysis.analysis import *

ROOT.TH1.SetDefaultSumw2(True)

f_iso = ROOT.TFile('output/important/output_2016_baseline_190814-iso-full.root')
hists_iso = TDirToNode(f_iso)
h_num_m = hists_iso['m']['baseline_m']['p0_pt']['W_J']
h_num_e = hists_iso['e']['baseline_e']['p0_pt']['W_J']
h_num_m.Add(h_num_e)

f_aiso = ROOT.TFile('output/important/output_2016_baseline_190814-iso-veto.root')
hists_aiso = TDirToNode(f_aiso)
h_den_m = hists_aiso['m']['baseline_m']['p0_pt']['W_J']
h_den_e = hists_aiso['e']['baseline_e']['p0_pt']['W_J']
h_den_m.Add(h_den_e)

h_num_m.Print("range")
h_den_m.Print("range")
h_num_m.Divide(h_den_m)
h_num_m.Print("range")

fout = ROOT.TFile('high_pt_photon_fakes_MC.root', 'RECREATE')



h2d = ROOT.TH2F('photon_fakes', '', 110, 100, 1200, 2, array('d', [0., 1.4442, 2.5]))
for i in xrange(1, h2d.GetNbinsX() + 1):
    x = h2d.GetXaxis().GetBinCenter(i)
    h2d.SetBinContent(i, 1, h_num_m.GetBinContent(h_num_m.GetXaxis().FindFixBin(x)))
    h2d.SetBinContent(i, 2, h_num_m.GetBinContent(h_num_m.GetXaxis().FindFixBin(x)))
    h2d.SetBinError(i, 1, h_num_m.GetBinError(h_num_m.GetXaxis().FindFixBin(x)))
    h2d.SetBinError(i, 2, h_num_m.GetBinError(h_num_m.GetXaxis().FindFixBin(x)))
h2d.Print('range')

ROOT.gDirectory.WriteTObject(h2d, 'photon_fakes')

fout.Close()