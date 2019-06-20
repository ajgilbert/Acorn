from ROOT import *

th1file =TFile.Open("histos_antiisosideband_mm17_0806.root")
antiiso_sideband = th1file.Get("data_H_mass_onrho_iso")

signal_file = TFile.Open("hzrho_sig_2017_1206_finest.root")
signal_hist = signal_file.Get("hzrho_mm")

w = RooWorkspace("w")


#x = RooRealVar("x","x",118,170)
#p1 = RooRealVar('p1', 'p1', 6.845, -10.0, 5000.0)
#p2 = RooRealVar('p2', 'p2', 6.7807, -10.0, 1000.0)
#p3 = RooRealVar('p3', 'p3', 0.30685, -10.0, 200.0)
#sqrtS = 13000.

#bkg = RooGenericPdf('bkg','(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))'%(sqrtS,sqrtS,sqrtS),RooArgList(x, p1, p2, p3))

#getattr(w,'import')(bkg)
w.factory("Chebychev::bkg(x[130,118,168],{c0[0.1,-1.,1.],c1[0.1,-1.,1.],c2[-0.1,-1.,1.]})")
#w.factory("Chebychev::bkg(x[130,118,168],{c0[-0.7349,-1.,1.],c1[0.108,-1.,1.],c2[0.0221,-1.,1.],c3[-0.01739,-1.,1.]})")
#w.factory("Exponential::bkg(x[130,118,168],coef[-0.030158,-1.,1.])")
w.factory("bkg_norm[11920.,1,100000]")
w.factory("SUM::addpdf(bkg_norm*bkg)")
#w.factory("SUM::bkg(bkg_norm[100,100000]*cheb)")
#w.factory("SUM::bkg(bkg_norm[100,100000]*expo)")
w.var("x").setRange("Range1",118,120)
w.var("x").setRange("Range2",130,168)

#w.factory("Bernstein::bkg(x[130,120,170],{coef1m17[0.152,0,1],coef2m17[0.059,0,1],coef3m17[0.043,0,1]})") rho window
data_obs=RooDataHist("data_obs","data_obs",RooArgList(w.var("x")),antiiso_sideband)
getattr(w,'import')(data_obs)
signal_dh=RooDataHist("signal_dh","signal_dh",RooArgList(w.var("x")),signal_hist)
#sig = RooHistPdf("sig","sig",RooArgSet(w.var("x")),signal_dh)
getattr(w,'import')(signal_dh)
#getattr(w,'import')(sig)
#w.factory("bkg_norm[40494.,1000,80000]")

pdf_bkg = w.pdf("addpdf")

dataset = w.data("data_obs")

fitres = pdf_bkg.fitTo(dataset,RooFit.SumW2Error(False),RooFit.Save(),RooFit.Range("Range1,Range2"))
fitres.Print()

#xframe = x.frame()
#data_obs.plotOn(xframe)
#pdf_bkg.plotOn(xframe,RooFit.Range("FullRange"))
#pdf_bkg.paramOn(xframe)
#bernstein.plotOn(xframe)

#extpdf.plotOn(xframe)
#extpdf.paramOn(xframe)
#c1=TCanvas("c1","c1")
#xframe.Draw()
#c1.SaveAs("wsfit.pdf")

w.writeToFile("workspace_mm_nominal.root")
w.Delete()

#th1file2 =TFile.Open("histos_antiisosideband_ee17.root")
#antiiso_sideband2 = th1file.Get("data_H_mass_onrho_antiiso2")
#antiiso_sideband2.Scale(19148./antiiso_sideband2.Integral())
#antiiso_sideband2.Scale(12393./antiiso_sideband2.Integral()) rho window

#signal_hist_ee = signal_file.Get("hzrho_ee")

#datahist = RooDataHist("datahist","datahist",RooArgList(x),antiiso_sideband)

#w = RooWorkspace("w")
#w.factory("Bernstein::bkg(x[130,120,170],{coef1e17[0.141,0,1],coef2e17[0.072,0,1],coef3e17[0.043,0,1]})")
#w.factory("Bernstein::bkg(x[130,120,170],{coef1e17[0.152,0,1],coef2e17[0.073,0,1],coef3e17[0.042,0,1]})") rho window
#data_obs=RooDataHist("data_obs","data_obs",RooArgList(w.var("x")),antiiso_sideband)
#getattr(w,'import')(data_obs)
#signal_dh=RooDataHist("signal_dh","signal_dh",RooArgList(w.var("x")),signal_hist_ee)
#sig = RooHistPdf("sig","sig",RooArgSet(w.var("x")),signal_dh)
#getattr(w,'import')(signal_dh)
#getattr(w,'import')(sig)
#w.factory("bkg_norm[19148.,1000,50000]")
#w.factory("bkg_norm[12393.,1000,50000]") rho window
#w.writeToFile("workspace_ee17.root")
#w.Delete()

