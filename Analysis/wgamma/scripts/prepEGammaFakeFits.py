import ROOT
from array import array
from Acorn.Analysis.analysis import *
import argparse

ROOT.RooWorkspace.imp = getattr(ROOT.RooWorkspace, 'import')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(0)

parser = argparse.ArgumentParser()
parser.add_argument('-y', '--year', default='2016')
parser.add_argument('-i', '--input', default='output_2016_electron_fakes.root')
parser.add_argument('--bkg-sub', type=float, default=0.0)
parser.add_argument('--label', default='nominal')
# parser.add_argument('-o', default='2016', choices=['2016', '2017', '2018'])
# parser.add_argument('--indir', default='output/130818/wgamma_2016_v2/WGamma/')
args = parser.parse_args()

sub_bkgs = ['DY_XZG_R', 'ZG_IZG_R', 'WG_R']

if args.bkg_sub:
    print '>> Subtracting backgrounds %s from data_obs' % str(sub_bkgs)

name = 'EGammaFakes'

fin = ROOT.TFile(args.input)

pt_bins = [30, 35, 40, 60, 100, 200]
eta_bins = [0, 1.4442, 2.5]

hist = ROOT.TH2F(name, name, len(pt_bins) - 1, array('d', pt_bins), len(eta_bins) - 1, array('d', eta_bins))
hist.GetXaxis().SetTitle('p_pt')
hist.GetYaxis().SetTitle('abs(p_eta)')

for sample in ['data_obs', 'DY_E']:
    outfile = ROOT.TFile('EFakesTP_%s_%s_%s.root' % (args.year, sample, args.label), 'RECREATE')
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


            if args.bkg_sub > 0.0 and sample == 'data_obs':
                h_pass_bkg = fin.Get('e/pass_%s/l0p0_M/%s' % (dirname, sub_bkgs[0]))
                h_fail_bkg = fin.Get('e/fail_%s/l0p0_M/%s' % (dirname, sub_bkgs[0]))
                for bkg in sub_bkgs[1:]:
                    h_pass_bkg.Add(fin.Get('e/pass_%s/l0p0_M/%s' % (dirname, bkg)))
                    h_fail_bkg.Add(fin.Get('e/fail_%s/l0p0_M/%s' % (dirname, bkg)))

                h_fail.Add(h_fail_bkg, -1. * args.bkg_sub)
                h_pass.Add(h_pass_bkg, -1. * args.bkg_sub)

                for ibin in xrange(1, h_fail.GetNbinsX() + 1):
                    if h_fail.GetBinContent(ibin) < 0.:
                        h_fail.SetBinContent(ibin, 0.)
                for ibin in xrange(1, h_pass.GetNbinsX() + 1):
                    if h_pass.GetBinContent(ibin) < 0.:
                        h_pass.SetBinContent(ibin, 0.)


            dat = wsp.imp(ROOT.RooDataHist(binname, '', ROOT.RooArgList(var),
                          ROOT.RooFit.Index(wsp.factory('cat[fail,pass]')),
                          ROOT.RooFit.Import('fail', h_fail),
                          ROOT.RooFit.Import('pass', h_pass))
                          )

    # rdh_pass_bkg = ROOT.RooDataHist('rdh_pass_bkg', '', ROOT.RooArgList(var), h_pass_bkg)
    # rdh_fail_bkg = ROOT.RooDataHist('rdh_fail_bkg', '', ROOT.RooArgList(var), h_fail_bkg)

    # wsp.imp(ROOT.RooHistPdf('backgroundPass', '', ROOT.RooArgSet(var), rdh_pass_bkg))
    # wsp.imp(ROOT.RooHistPdf('backgroundFail', '', ROOT.RooArgSet(var), rdh_fail_bkg))

    outfile.cd()
    wsp.Write()
    hist.Write()
    wsp.Delete()

    outfile.Close()
