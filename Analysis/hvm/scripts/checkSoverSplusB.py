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


#infile_data = ROOT.TFile.Open("output/PROD-29032018-17/SingleMuon.root")
#infile_mc = ROOT.TFile.Open("output/PROD-29032018-17/HZrhoMM.root")
infile_data = ROOT.TFile.Open("adddileptrig/PROD-29032018-17/SingleMuon.root")
infile_data_doub = ROOT.TFile.Open("adddileptrig/PROD-29032018-17/DoubleMuon.root")
infile_mc = ROOT.TFile.Open("adddileptrig/PROD-29032018-17/HZrhoMM_0.root")


intree_data = infile_data.Get("DiMuonMesonAnalysis")
intree_data_doub = infile_data_doub.Get("DiMuonMesonAnalysis")
intree_mc = infile_mc.Get("DiMuonMesonAnalysis")


data_H_mass_antiiso = ROOT.TH1F("data_H_mass_antiiso","data_H_mass_antiiso",20,120,140)
data_H_mass_antiiso.Sumw2()
mc_H_mass = ROOT.TH1F("mc_H_mass","mc_H_mass",20,120,140)
mc_H_mass.Sumw2()


#intree_data.Draw("highestpt_pair_reco_higgs_mass>>data_H_mass_antiiso","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso>1.0&&highestpt_pair_looser_iso<2.0&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&(highestpt_pair_mass>0.5&&highestpt_pair_mass<1.)")
#intree_mc.Draw("highestpt_pair_reco_higgs_mass>>mc_H_mass","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&(highestpt_pair_mass>0.5&&highestpt_pair_mass<1.)")
intree_data.Draw("highestpt_pair_reco_higgs_mass>>data_H_mass_antiiso","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso>1.0&&highestpt_pair_looser_iso<2.0&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1>30")
#intree_data_doub.Draw("highestpt_pair_reco_higgs_mass>>+data_H_mass_antiiso","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso>1.0&&highestpt_pair_looser_iso<2.0&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1<30")
intree_mc.Draw("highestpt_pair_reco_higgs_mass>>mc_H_mass","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)")


print "data yield in anti-isolated region: ", data_H_mass_antiiso.Integral()
#print "data yield in isolated region: ",intree_data.GetEntries("highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&(highestpt_pair_mass>0.5&&highestpt_pair_mass<1.)")
print "data yield in isolated region (SL): ",intree_data.GetEntries("highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1>30")
print "data yield in isolated region (DL): ",intree_data_doub.GetEntries("highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1<30")
#data_H_mass_antiiso.Scale(31915./data_H_mass_antiiso.Integral())
data_H_mass_antiiso.Scale(31450./data_H_mass_antiiso.Integral())
#data_H_mass_antiiso.Scale(31965./data_H_mass_antiiso.Integral())
#data_H_mass_antiiso.Scale(14648.0/data_H_mass_antiiso.Integral())


ams_binned=0;
for i in range(1,21):
  tmp_ams=ROOT.RooStats.AsimovSignificance(mc_H_mass.GetBinContent(i),data_H_mass_antiiso.GetBinContent(i))
  ams_binned+=tmp_ams*tmp_ams

print "sqrt of summed sqrt binned ams: ", ROOT.TMath.Sqrt(ams_binned)

mc_H_mass.SetMarkerColor(ROOT.kGreen+3)
mc_H_mass.Scale(41519./100000.)
canvas = ROOT.TCanvas("c1","c1")
pads = plot.OnePad()
axish = createAxisHists(1,data_H_mass_antiiso,data_H_mass_antiiso.GetXaxis().GetXmin(),data_H_mass_antiiso.GetXaxis().GetXmax()-0.01)
#pads[0].SetLogy()
axish[0].GetYaxis().SetRangeUser(1,2*data_H_mass_antiiso.GetMaximum())
axish[0].GetXaxis().SetTitle("m_{#mu#mu#rho} [GeV]")
axish[0].GetYaxis().SetTitle("Events")
pads[0].cd()
axish[0].Draw()
data_H_mass_antiiso.Draw("PSAME")
mc_H_mass.Draw("PSAME")
legend = plot.PositionedLegend(0.4, 0.1, 3, 0.015)
legend.AddEntry(data_H_mass_antiiso, "Background", "P")
legend.AddEntry(mc_H_mass, "H#rightarrowZ(#mu#mu)#rho (1pb)", "P")
legend.Draw("SAME")
axish[0].Draw("axissame")
plot.DrawTitle(pads[0], "41.5 fb^{-1} (13 TeV)", 3)
canvas.SaveAs("mH_antiiso_Z60to120_MCiso_mm_trigcocktail.pdf")

infile_data_ee = ROOT.TFile.Open("adddileptrig/PROD-29032018-17/SingleElectron.root")
infile_data_doub_ee = ROOT.TFile.Open("adddileptrig/PROD-29032018-17/DoubleEG.root")
infile_mc_ee = ROOT.TFile.Open("adddileptrig/PROD-29032018-17/HZrhoEE_0.root")

intree_data_ee = infile_data_ee.Get("DiElectronMesonAnalysis")
intree_data_doub_ee = infile_data_doub_ee.Get("DiElectronMesonAnalysis")
intree_mc_ee = infile_mc_ee.Get("DiElectronMesonAnalysis")


data_H_mass_antiiso_ee = ROOT.TH1F("data_H_mass_antiiso_ee","data_H_mass_antiiso_ee",20,120,140)
data_H_mass_antiiso_ee.Sumw2()
mc_H_mass_ee = ROOT.TH1F("mc_H_mass_ee","mc_H_mass_ee",20,120,140)
mc_H_mass_ee.Sumw2()


#intree_data.Draw("highestpt_pair_reco_higgs_mass>>data_H_mass_antiiso","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso>1.0&&highestpt_pair_looser_iso<2.0&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&(highestpt_pair_mass>0.5&&highestpt_pair_mass<1.)")
#intree_mc.Draw("highestpt_pair_reco_higgs_mass>>mc_H_mass","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&(highestpt_pair_mass>0.5&&highestpt_pair_mass<1.)")
intree_data_ee.Draw("highestpt_pair_reco_higgs_mass>>data_H_mass_antiiso_ee","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso>1.0&&highestpt_pair_looser_iso<2.0&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1>38")
#intree_data_doub_ee.Draw("highestpt_pair_reco_higgs_mass>>+data_H_mass_antiiso_ee","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso>1.0&&highestpt_pair_looser_iso<2.0&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1<38")
intree_mc_ee.Draw("highestpt_pair_reco_higgs_mass>>mc_H_mass_ee","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)")


print "data yield in anti-isolated region: ", data_H_mass_antiiso_ee.Integral()
#print "data yield in isolated region: ",intree_data.GetEntries("highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&(highestpt_pair_mass>0.5&&highestpt_pair_mass<1.)")
print "data yield in isolated region (SL): ",intree_data_ee.GetEntries("highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1>38")
print "data yield in isolated region (DL): ",intree_data_doub_ee.GetEntries("highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_1_pt>10||highestpt_pair_2_pt>10)&&pt_1<38")

data_H_mass_antiiso_ee.Scale(14159./data_H_mass_antiiso_ee.Integral())
#data_H_mass_antiiso_ee.Scale(15140./data_H_mass_antiiso_ee.Integral())
#data_H_mass_antiiso_ee.Scale(14952./data_H_mass_antiiso_ee.Integral())
#data_H_mass_antiiso.Scale(14648.0/data_H_mass_antiiso.Integral())

ams_binned_ee=0;
for i in range(1,21):
  tmp_ams=ROOT.RooStats.AsimovSignificance(mc_H_mass_ee.GetBinContent(i),data_H_mass_antiiso_ee.GetBinContent(i))
  ams_binned+=tmp_ams*tmp_ams

print "sqrt of summed sqrt binned ams, ee channel: ", ROOT.TMath.Sqrt(ams_binned)

mc_H_mass_ee.SetMarkerColor(ROOT.kGreen+3)
mc_H_mass_ee.Scale(41519./150000.)
canvas1 = ROOT.TCanvas("c2","c2")
pads1 = plot.OnePad()
axish1 = createAxisHists(1,data_H_mass_antiiso_ee,data_H_mass_antiiso_ee.GetXaxis().GetXmin(),data_H_mass_antiiso_ee.GetXaxis().GetXmax()-0.01)
#pads[0].SetLogy()
axish1[0].GetYaxis().SetRangeUser(1,2*data_H_mass_antiiso_ee.GetMaximum())
axish1[0].GetXaxis().SetTitle("m_{ee#rho} [GeV]")
axish1[0].GetYaxis().SetTitle("Events")
pads1[0].cd()
axish1[0].Draw()
data_H_mass_antiiso_ee.Draw("PSAME")
mc_H_mass_ee.Draw("PSAME")
legend1 = plot.PositionedLegend(0.4, 0.1, 3, 0.015)
legend1.AddEntry(data_H_mass_antiiso_ee, "Background", "P")
legend1.AddEntry(mc_H_mass_ee, "H#rightarrow Z(ee)#rho (1pb)", "P")
legend1.Draw("SAME")
axish1[0].Draw("axissame")
plot.DrawTitle(pads1[0], "41.5 fb^{-1} (13 TeV)", 3)
canvas1.SaveAs("mH_antiiso_Z60to120_MCiso_ee_trigcocktail.pdf")
