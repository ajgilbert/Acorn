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


infile_data = ROOT.TFile.Open("output-2505/PROD-03052019-18/SingleMuon.root")
infile_DYMC = ROOT.TFile.Open("output-2505/PROD-03052019-18/DYJetsToLL_M-50-madgraphMLM.root")
intree_data = infile_data.Get("DiMuonMesonAnalysis")
intree_DYMC = infile_DYMC.Get("DiMuonMesonAnalysis")


data_H_mass_onrho_iso = ROOT.TH1F("data_H_mass_onrho_iso","data_H_mass_onrho_iso",70,110,180)
data_H_mass_onrho_iso.Sumw2()
intree_data.Draw("highestpt_pair_reco_higgs_mass>>data_H_mass_onrho_iso","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&(highestpt_pair_reco_higgs_mass<120||highestpt_pair_reco_higgs_mass>130)&&highestpt_pair_mass>0.5&&highestpt_pair_mass<1.0")
data_H_mass_onrho_iso.SetMarkerStyle(20)
data_H_mass_onrho_iso.SetMarkerColor(ROOT.kBlack)

data_H_mass_onrho_antiiso2 = ROOT.TH1F("data_H_mass_onrho_antiiso2","data_H_mass_onrho_antiiso2",70,110,180)
data_H_mass_onrho_antiiso2.Sumw2()
intree_data.Draw("highestpt_pair_reco_higgs_mass>>data_H_mass_onrho_antiiso2","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<2.&&highestpt_pair_looser_iso>1.&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&highestpt_pair_mass>0.5&&highestpt_pair_mass<1.0")
data_H_mass_onrho_antiiso2.SetMarkerStyle(20)
data_H_mass_onrho_antiiso2.SetMarkerColor(ROOT.kBlack)
#print "Nentries in SR 120-140:"
#print intree_data.GetEntries("highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&(trg_1||trg_2)&&m_ll>80&&m_ll<100&&highestpt_pair_reco_higgs_mass>120&&highestpt_pair_reco_higgs_mass<140")

dy_H_mass_onrho_iso = ROOT.TH1F("dy_H_mass_onrho_iso","dy_H_mass_onrho_iso",70,110,180)
dy_H_mass_onrho_iso.Sumw2()
intree_DYMC.Draw("highestpt_pair_reco_higgs_mass>>dy_H_mass_onrho_iso","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso<0.5&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&highestpt_pair_mass>0.5&&highestpt_pair_mass<1.0")
dy_H_mass_onrho_iso.SetMarkerStyle(20)
dy_H_mass_onrho_iso.SetMarkerColor(ROOT.kBlack)


dy_H_mass_onrho_antiiso2 = ROOT.TH1F("dy_H_mass_onrho_antiiso2","dy_H_mass_onrho_antiiso2",70,110,180)
dy_H_mass_onrho_antiiso2.Sumw2()
intree_DYMC.Draw("highestpt_pair_reco_higgs_mass>>dy_H_mass_onrho_antiiso2","highestpt_pair_looser_iso>-1&&highestpt_pair_looser_iso>1.&&highestpt_pair_looser_iso<2.0&&(trg_1||trg_2)&&m_ll>60&&m_ll<120&&highestpt_pair_mass>0.5&&highestpt_pair_mass<1.0")
dy_H_mass_onrho_antiiso2.SetMarkerStyle(20)
dy_H_mass_onrho_antiiso2.SetMarkerColor(ROOT.kBlack)

outfile = ROOT.TFile.Open("histos_antiisosideband_mumu18_rhocut.root","RECREATE")
data_H_mass_onrho_iso.Write()
data_H_mass_onrho_antiiso2.Write()
dy_H_mass_onrho_antiiso2.Write()
dy_H_mass_onrho_iso.Write()
outfile.Close()
