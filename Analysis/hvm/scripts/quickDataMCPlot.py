import ROOT
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
#from Acorn.Analysis.analysis import *
from array import array
import argparse

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
#ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

def createAxisHists(n,src,xmin=0,xmax=499):
  result = []
  for i in range(0,n):
    res = src.Clone()
    res.Reset()
    res.SetTitle("")
    res.SetName("axis%(i)d"%vars())
    res.SetAxisRange(xmin,xmax)
    res.SetStats(0)
    result.append(res)
  return result


infile_data = ROOT.TFile.Open("output/PROD-29032018-17/SingleElectron.root")
#infile_data = ROOT.TFile.Open("output/PROD-29032018-17/SingleMuon.root")
infile_DYMC = ROOT.TFile.Open("output/PROD-29032018-17/DYJetsToLL_M-50-amcNLO.root")
#infile_DYMC = ROOT.TFile.Open("output/PROD-29032018-17/DYJetsToLL_M-50-madgraphMLM.root")

intree_data = infile_data.Get("DiElectronMesonAnalysis")
intree_DYMC = infile_DYMC.Get("DiElectronMesonAnalysis")
#intree_DYMCNLO = infile_DYMCNLO.Get("DiMuonMesonAnalysis")
LOMC_counter = infile_DYMC.Get("counters")

data_Z_mass = ROOT.TH1F("data_Z_mass","data_Z_mass",60,60,120)
data_Z_mass.Sumw2()
intree_data.Draw("m_ll>>data_Z_mass","(trg_1||trg_2)&&dr_ll>0.5")
data_Z_mass.SetMarkerStyle(20)
data_Z_mass.SetMarkerColor(ROOT.kBlack)

data_Z_pt = ROOT.TH1F("data_Z_pt","data_Z_pt",60,0,120)
data_Z_pt.Sumw2()
intree_data.Draw("pt_ll>>data_Z_pt","(trg_1||trg_2)&&dr_ll>0.5")
data_Z_pt.SetMarkerStyle(20)
data_Z_pt.SetMarkerColor(ROOT.kBlack)


LOMC_Z_mass = ROOT.TH1F("LOMC_Z_mass","LOMC_Z_mass",60,60,120)
LOMC_Z_mass.SetMarkerSize(0)
LOMC_Z_mass.Sumw2()
#intree_DYMC.Draw("pt_1>>LOMC_Z_mass","(trg_1||trg_2)*(wt_pu)")
#intree_DYMC.Draw("m_ll>>LOMC_Z_mass","((trg_1||trg_2)&&dr_ll>0.5)*(wt_pu)")
intree_DYMC.Draw("m_ll>>LOMC_Z_mass","((trg_1||trg_2)&&dr_ll>0.5)*(wt_pu*wt_1*wt_2*wt_trg*gen_wt)")
LOMC_Z_mass.Scale(6225.42*41519.180032466/LOMC_counter.GetBinContent(2))
LOMC_Z_mass.SetFillStyle(1001)
LOMC_Z_mass.SetFillColor(ROOT.TColor.GetColor(100,192,232))
LOMC_Z_mass.SetLineColor(ROOT.kBlack)
LOMC_Z_mass_unc = LOMC_Z_mass.Clone("LOMC_Z_mass_unc")
LOMC_Z_mass_unc.SetFillColor(plot.CreateTransparentColor(12,0.4))
LOMC_Z_mass_unc.SetLineColor(0)
LOMC_Z_mass_unc.SetMarkerSize(0)

LOMC_Z_pt = ROOT.TH1F("LOMC_Z_pt","LOMC_Z_pt",60,0,120)
LOMC_Z_pt.SetMarkerSize(0)
LOMC_Z_pt.Sumw2()
#intree_DYMC.Draw("pt_2>>LOMC_Z_pt","(trg_1||trg_2)*(wt_pu)")
#intree_DYMC.Draw("pt_ll>>LOMC_Z_pt","((trg_1||trg_2)&&dr_ll>0.5)*(wt_pu)")
intree_DYMC.Draw("pt_ll>>LOMC_Z_pt","((trg_1||trg_2)&&dr_ll>0.5)*(wt_pu*wt_1*wt_2*wt_trg*gen_wt)")
LOMC_Z_pt.Scale(6225.42*41519.180032466/LOMC_counter.GetBinContent(2))
LOMC_Z_pt.SetFillStyle(1001)
LOMC_Z_pt.SetFillColor(ROOT.TColor.GetColor(100,192,232))
LOMC_Z_pt.SetLineColor(ROOT.kBlack)
LOMC_Z_pt_unc = LOMC_Z_pt.Clone("LOMC_Z_pt_unc")
LOMC_Z_pt_unc.SetFillColor(plot.CreateTransparentColor(12,0.4))
LOMC_Z_pt_unc.SetLineColor(0)
LOMC_Z_pt_unc.SetMarkerSize(0)


canvas = ROOT.TCanvas("c1","c1")
pads = plot.TwoPadSplit(0.29,0.01,0.01)
axish = createAxisHists(2,LOMC_Z_mass,LOMC_Z_mass.GetXaxis().GetXmin(),LOMC_Z_mass.GetXaxis().GetXmax()-0.01)
pads[0].SetLogy()
axish[0].GetYaxis().SetRangeUser(1,12*LOMC_Z_mass.GetMaximum())
axish[0].GetXaxis().SetTitleSize(0)
axish[0].GetXaxis().SetLabelSize(0)
axish[1].GetYaxis().SetRangeUser(0.8,1.2)
axish[1].GetXaxis().SetTitle("m_{ee} [GeV]")
axish[0].GetYaxis().SetTitle("Events")
axish[1].GetYaxis().SetTitle("Obs/Exp")
ratio_mchist = plot.MakeRatioHist(LOMC_Z_mass_unc, LOMC_Z_mass_unc, True, False)
ratio_datahist = plot.MakeRatioHist(data_Z_mass, LOMC_Z_mass, True, False)

pads[0].cd()
axish[0].Draw()
LOMC_Z_mass.Draw("HISTSAME")
LOMC_Z_mass_unc.Draw("e2same")
data_Z_mass.Draw("PSAME")
legend = plot.PositionedLegend(0.2, 0.1, 3, 0.015)
legend.AddEntry(LOMC_Z_mass, "Z#rightarrow ee","F")
legend.AddEntry(data_Z_mass, "Data", "P")
legend.Draw("SAME")
axish[0].Draw("axissame")
pads[1].cd()
axish[1].Draw("axis")
ratio_mchist.SetMarkerSize(0)
ratio_mchist.Draw("e2same")
ratio_datahist.Draw("e0same")
pads[1].RedrawAxis("G")
plot.DrawTitle(pads[0], "41.5 fb^{-1} (13 TeV-2017)", 3)

canvas.SaveAs("mee_inclusive_nlo_weight_2017.pdf")
#canvas.SaveAs("mmumu_inclusive_weight.pdf")

canvas2 = ROOT.TCanvas("c1","c1")
pads2 = plot.TwoPadSplit(0.29,0.01,0.01)
axish2 = createAxisHists(2,LOMC_Z_pt,LOMC_Z_pt.GetXaxis().GetXmin(),LOMC_Z_pt.GetXaxis().GetXmax()-0.01)
pads2[0].SetLogy()
axish2[0].GetYaxis().SetRangeUser(1,12*LOMC_Z_pt.GetMaximum())
axish2[0].GetXaxis().SetTitleSize(0)
axish2[0].GetXaxis().SetLabelSize(0)
axish2[1].GetYaxis().SetRangeUser(0.8,1.2)
axish2[1].GetXaxis().SetTitle("p_{T}(ee) [GeV]")
axish2[0].GetYaxis().SetTitle("Events")
axish2[1].GetYaxis().SetTitle("Obs/Exp")
ratio_mchist_pt = plot.MakeRatioHist(LOMC_Z_pt_unc, LOMC_Z_pt_unc, True, False)
ratio_datahist_pt = plot.MakeRatioHist(data_Z_pt, LOMC_Z_pt, True, False)

pads2[0].cd()
axish2[0].Draw()
LOMC_Z_pt.Draw("HISTSAME")
LOMC_Z_pt_unc.Draw("e2same")
data_Z_pt.Draw("PSAME")
legend2 = plot.PositionedLegend(0.2, 0.1, 3, 0.015)
legend2.AddEntry(LOMC_Z_pt, "Z#rightarrowee","F")
legend2.AddEntry(data_Z_pt, "Data", "P")
legend2.Draw("SAME")
axish2[0].Draw("axissame")
pads2[1].cd()
axish2[1].Draw("axis")
ratio_mchist_pt.SetMarkerSize(0)
ratio_mchist_pt.Draw("e2same")
ratio_datahist_pt.Draw("e0same")
pads2[1].RedrawAxis("G")
plot.DrawTitle(pads2[0], "41.5 fb^{-1} (13 TeV-2017)", 3)

canvas2.SaveAs("ptee_inclusive_nlo_weight_2017.pdf")
#canvas2.SaveAs("ptmumu_inclusive_weight.pdf")

