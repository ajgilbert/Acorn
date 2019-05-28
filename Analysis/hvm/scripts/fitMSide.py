from ROOT import *
gSystem.Load("libHiggsAnalysisCombinedLimit.so")

#th1file = TFile.Open("histos_antiisosideband_ee17_rhocut.root")
th1file = TFile.Open("histos_antiisosideband_mumu17_rhocut.root")
antiiso_sideband = th1file.Get("data_H_mass_onrho_iso")
#print antiiso_sideband.Integral()
#antiiso_sideband.Scale(11920./antiiso_sideband.Integral())

x = RooRealVar("x","x",118,170)
x.setRange("Range1",118,120)
x.setRange("Range2",130,170)
x.setRange("FullRange",118,170)

datahist = RooDataHist("datahist","datahist",RooArgList(x),antiiso_sideband)

coef1 = RooRealVar("coef1","coef1",0.5,-10,10)
coef2 = RooRealVar("coef2","coef2",0.5,-10,10)
coef3 = RooRealVar("coef3","coef3",0.5,-10,10)

c0 = RooRealVar("c0","coefficient #0",0.1,-1.,1)
c1 = RooRealVar("c1","coefficient #1",0.1,-1.,1)
c2 = RooRealVar("c2","coefficient #2",-0.1,-1.,1)
c3 = RooRealVar("c3","coefficient #2",-0.1,-1.,1)

cl0 = RooRealVar("cl0","cl0",-0.1,-1.,1.)
cl1 = RooRealVar("cl1","cl1",-0.1,-1.,1.)

coefb1 = RooRealVar("coefb1","coefb1",0.1,0,1)
coefb2 = RooRealVar("coefb2","coefb2",0.1,0,1)
coefb3 = RooRealVar("coefb3","coefb3",0.1,0,1)
coefb4 = RooRealVar("coefb4","coefb4",0.1,0,1)
coefb5 = RooRealVar("coefb5","coefb5",0.1,0,1)
coefb6 = RooRealVar("coefb6","coefb6",0.1,0,1)
coefb7 = RooRealVar("coefb7","coefb7",0.1,0,1)

exp1 = RooRealVar("exp1","exp1",-5,-10,-0.1)
exp2 = RooRealVar("exp2","exp2",-3,-10,-0.1)
exp3 = RooRealVar("exp2","exp2",-3,-10,-0.1)

l1 = RooRealVar("l1","l1",-4)
l2 = RooRealVar("l2","l2",-3)
l3 = RooRealVar("l3","l3",-5)


bernstein = RooBernstein("bernstein","bernstein",x,RooArgList(coefb1,coefb2,coefb3))
cheby = RooChebychev("chebychev","chebychev",x,RooArgList(c0))
chebyprod = RooPolynomial("chebyprod","chebyrod",x,RooArgList(c0,c1))
expo1 = RooExponential("expo1","expo1",x,coef1)
expo2 = RooExponential("expo2","expo2",x,coef2)
expo3 = RooExponential("expo3","expo3",x,coef3)
expoprod = RooProdPdf("expoprod","expoprod",RooArgList(chebyprod,expo1))
pow1 = RooPower("pow1","pow1",x,exp1)
pow2 = RooPower("pow2","pow2",x,exp2)
lau1 = RooPower("lau1","lau1",x,l1)
lau2 = RooPower("lau2","lau2",x,l2)
lau3 = RooPower("lau3","lau3",x,l3)
normterm = RooRealVar("normterm","normterm",11920,1.,100000)
normterm1 = RooRealVar("normterm1","normterm1",11920,1.,100000)
normterm2 = RooRealVar("normterm2","normterm2",11920,1.,100000)
normterm3 = RooRealVar("normterm3","normterm3",11920,1.,100000)
#addpdf =RooExtendPdf("addpdf","addpdf",bernstein,normterm)
#addpdf = RooAddPdf("addpdf","addpdf",RooArgList(bernstein),RooArgList(normterm))
#addpdf = RooAddPdf("addpdf","addpdf",RooArgList(expoprod),RooArgList(normterm))
addpdf = RooAddPdf("addpdf","addpdf",RooArgList(expo1,expo2,expo3),RooArgList(normterm1,normterm2,normterm3))
#addpdf = RooAddPdf("addpdf","addpdf",RooArgList(pow1),RooArgList(normterm1))
#addpdf = RooAddPdf("addpdf","addpdf",RooArgList(lau1),RooArgList(normterm1))
#addpdf = RooAddPdf("addpdf","addpdf",RooArgList(expo1,expo2),RooArgList(normterm1,normterm2))

#fitresult = extpdf.fitTo(datahist,RooFit.Extended(True),RooFit.SumW2Error(False),RooFit.Save())
#fitresult = addpdf.fitTo(datahist,RooFit.Extended(True),RooFit.SumW2Error(False),RooFit.Save())
fitresult = addpdf.fitTo(datahist,RooFit.SumW2Error(False),RooFit.Save(),RooFit.Range("Range1,Range2"))
fitresult.Print()
#print fitresult.minNll()

print datahist.numEntries()
rss = 0
for i in range(0,datahist.numEntries()):
  if i < 2 or i > 11:
    x.setVal(antiiso_sideband.GetBinCenter(9+i))
    #rss += pow((normterm.getVal()*addpdf.getVal(RooArgSet(x))-datahist.weight(datahist.get(i))),2)
    rss += pow((normterm1.getVal()*pow1.getVal(RooArgSet(x))+normterm2.getVal()*pow2.getVal(RooArgSet(x))-datahist.weight(datahist.get(i))),2)
  #rss += pow((1.8373e+04*expo1.getVal(RooArgSet(x))+1.8431e+04*expo2.getVal(RooArgSet(x))-datahist.weight(datahist.get(i))),2)
print "RSS: ", rss

xframe = x.frame()
coefb1.Print()
coefb2.Print()
coefb3.Print()
normterm.getVal(RooArgSet(x))
datahist.plotOn(xframe)
normterm.getVal(RooArgSet(x))
coefb1.Print()
coefb2.Print()
coefb3.Print()
#addpdf.forceNumInt()
addpdf.plotOn(xframe,RooFit.Range("FullRange"))
#bernstein.plotOn(xframe)

addpdf.paramOn(xframe)
#extpdf.plotOn(xframe)
#extpdf.paramOn(xframe)
c1=TCanvas("c1","c1")
xframe.Draw()
print "Chi2/ndof: ",xframe.chiSquare(fitresult.floatParsFinal().getSize())
print "Number of floating parameters: ",fitresult.floatParsFinal().getSize()
c1.SaveAs("Fit2017_iso_rhocut_mm_expoprod.pdf")

