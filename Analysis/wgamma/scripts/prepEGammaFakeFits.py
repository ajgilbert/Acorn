import ROOT
import glob
import sys
# import json
from array import array
from Acorn.Analysis.analysis import *

ROOT.RooWorkspace.imp = getattr(ROOT.RooWorkspace, 'import')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(0)

name = 'EGammaFakes'

fin = ROOT.TFile('output_2016_electron_fakes.root')

pt_bins = [30, 35, 40, 60, 100, 200]
eta_bins = [0, 1.4442, 2.5]

hist = ROOT.TH2F(name, name, len(pt_bins) - 1, array('d', pt_bins), len(eta_bins) - 1, array('d', eta_bins))
hist.GetXaxis().SetTitle('p_pt')
hist.GetYaxis().SetTitle('abs(p_eta)')

for sample in ['data_obs', 'DY_E']:
    outfile = ROOT.TFile('ZeeTP_%s.root' % sample, 'RECREATE')
    wsp = ROOT.RooWorkspace('wsp_'+name, '')
    var = wsp.factory('m_ll[100,75,125]')


    for ieta in range(1, hist.GetNbinsY() + 1):
        for ipt in range(1, hist.GetNbinsX() + 1):
            binname = '%s>=%g && %s<%g && %s>=%g && %s<%g' % (
                    'p_pt', hist.GetXaxis().GetBinLowEdge(ipt),
                    'p_pt', hist.GetXaxis().GetBinUpEdge(ipt),
                    'abs(p_eta)', hist.GetYaxis().GetBinLowEdge(ieta),
                    'abs(p_eta)', hist.GetYaxis().GetBinUpEdge(ieta),
                    )
            dirname = '%g_%g_%g_%g' % (hist.GetYaxis().GetBinLowEdge(ieta), hist.GetYaxis().GetBinUpEdge(ieta), hist.GetXaxis().GetBinLowEdge(ipt), hist.GetXaxis().GetBinUpEdge(ipt))
            dirname = dirname.replace('.', 'p')
            print dirname
            h_fail = fin.Get('e/fail_%s/l0p0_M/%s' % (dirname, sample))
            h_pass = fin.Get('e/pass_%s/l0p0_M/%s' % (dirname, sample))

            h_pass_bkg = fin.Get('e/pass_%s/l0p0_M/WG' % dirname)
            h_pass_bkg.Add(fin.Get('e/pass_%s/l0p0_M/DY_R' % dirname))
            h_fail_bkg = fin.Get('e/fail_%s/l0p0_M/WG' % dirname)
            h_fail_bkg.Add(fin.Get('e/fail_%s/l0p0_M/DY_R' % dirname))

            if sample == 'data_obs':
                h_fail.Add(h_fail_bkg, -1)
                h_pass.Add(h_pass_bkg, -1)
    #     for b in cfg['bins']:
    #         hists[cfg['name']]['%s:fail' % b] = Hist('TH1F', sample=sample, var=[cfg['var']], binning=cfg['binning'], sel='((%s) && !(%s) && (%s))' % (b, cfg['probe'], cfg['tag']), wt='wt_def')
    #         hists[cfg['name']]['%s:pass' % b] = Hist('TH1F', sample=sample, var=[cfg['var']], binning=cfg['binning'], sel='((%s) && (%s) && (%s))' % (b, cfg['probe'], cfg['tag']), wt='wt_def')
    #         # drawlist.append((cfg['var'], '((%s) && !(%s) && (%s)) * wt' % (b, cfg['probe'], cfg['tag'])))
    #         # drawlist.append((cfg['var'], '((%s) && (%s) && (%s)) * wt' % (b, cfg['probe'], cfg['tag'])))

    # MultiDraw(hists, samples, 'WGTagAndProbe', mt_cores=4)

    # # hists = trees[sample].Draw(drawlist, compiled=True)
   

            dat = wsp.imp(ROOT.RooDataHist(binname, '', ROOT.RooArgList(var),
                          ROOT.RooFit.Index(wsp.factory('cat[fail,pass]')),
                          ROOT.RooFit.Import('fail', h_fail),
                          ROOT.RooFit.Import('pass', h_pass))
                          )

    # rdh_pass_bkg = ROOT.RooDataHist('rdh_pass_bkg', '', ROOT.RooArgList(var), h_pass_bkg)
    # rdh_fail_bkg = ROOT.RooDataHist('rdh_fail_bkg', '', ROOT.RooArgList(var), h_fail_bkg)

    # wsp.imp(ROOT.RooHistPdf('backgroundPass', '', ROOT.RooArgSet(var), rdh_pass_bkg))
    # print h_pass_bkg.Integral()
    # wsp.imp(ROOT.RooHistPdf('backgroundFail', '', ROOT.RooArgSet(var), rdh_fail_bkg))

    outfile.cd()
    wsp.Write()
    hist.Write()
    wsp.Delete()

    outfile.Close()
