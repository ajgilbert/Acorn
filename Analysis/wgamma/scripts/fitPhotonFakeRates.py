import ROOT
import CombineHarvester.CombineTools.plotting as plot
import math
import argparse
from array import array
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

ROOT.TH1.SetDefaultSumw2(True)
plot.ModTDRStyle(height=300)

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', default='')
parser.add_argument('-o', '--output', default='output_2016_photon_fakes_ratios_m.root')
parser.add_argument('--year', default='2016')
parser.add_argument('--var', default='p0_pt', choices=['p0_pt', 'p0_chiso'])

args = parser.parse_args()


def RebinHist(hist, bins):
    x = hist.Rebin(len(bins) - 1, "", array('d', bins))
    x.Copy(hist)


def To1Bin(hist):
    xmin = hist.GetXaxis().GetBinLowEdge(1)
    xmax = hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX())
    newhist = hist.Clone()
    RebinHist(newhist, [xmin, xmax])
    return newhist


def IntegralRatio(h_num, h_den):
    hr_num = To1Bin(h_num)
    hr_den = To1Bin(h_den)
    hr_num.Divide(hr_den)
    return (hr_num.GetBinContent(1), hr_num.GetBinError(1))


def Chi2Compatibility(*measurements):
    num = 0
    den = 0
    for x in measurements:
        num += (x[0] / (x[1] * x[1]))
        den += (1.0 / (x[1] * x[1]))
    weighted_mean = num / den
    standard_dev = math.sqrt(1.0 / den)
    # print (weighted_mean, standard_dev)
    chi2 = 0
    for x in measurements:
        chi2 += (math.pow(weighted_mean - x[0], 2) / math.pow(x[1], 2))
    pval = ROOT.Math.chisquared_cdf_c(chi2, len(measurements) - 1)
    return pval, weighted_mean, standard_dev

fin = ROOT.TFile(args.input)
hists = TDirToNode(fin)

var = args.var

channels = ['m', 'e']
selections = ['sig_l', 'sig_t', 'iso_t_sig_t', 'iso_t_sig_l']

common_syst = 0.1  # 10% to cover non-closure

config = {
    'barrel1': [
        ['30_40', ['30_40'], 0.85],
        ['40_50', ['40_50']],
        ['50_70', ['50_70']],
        ['70_100', ['70_100']],
        ['100_150', ['100_150']],
        ['150_300', ['150_300']]
    ],
    'barrel2': [
        ['30_40', ['30_40']],
        ['40_50', ['40_50']],
        ['50_70', ['50_70']],
        ['70_100', ['70_100']],
        ['100_300', ['100_150', '150_300']]
    ],
    'endcap1': [
        ['30_40', ['30_40']],
        ['40_70', ['40_50', '50_70']],
        ['70_300', ['70_100', '100_150', '150_300']],
    ],
    'endcap2': [
        ['30_40', ['30_40']],
        ['40_300', ['40_50', '50_70','70_100', '100_150', '150_300']],
    ],
    'barrel': [
        ['30_40', ['30_40']],
        ['40_50', ['40_50']],
        ['50_70', ['50_70']],
        ['70_100', ['70_100']],
        ['100_150', ['100_150']],
        ['150_300', ['150_300']]
    ],
    'endcap': [
        ['30_40', ['30_40']],
        ['40_70', ['40_50', '50_70']],
        ['70_300', ['70_100', '100_150', '150_300']],
    ],
}

eta_ranges = {
    'barrel1': ('0', '1.0'),
    'barrel2': ('1.0', '1.4442'),
    'endcap1': ('1.556', '2.1'),
    'endcap2': ('2.1', '2.5'),
}


reslines = []

eta_regions = ['barrel1', 'barrel2', 'endcap1', 'endcap2']

allresults = list()
bkg_subtract = [1.0, 0.8, 1.2]

ROOT.gStyle.SetOptFit(0)

for bkg_sub in bkg_subtract:
    allresults.append(list())
    all_fit_over_true = []

    for region in eta_regions:
        allresults[-1].append(list())
        for pt_range in config[region]:
            post = '_pt_%s' % pt_range[0]
            for sel in selections:
                seldir_cmb = '%s_%s%s' % (region, sel, post)
                print seldir_cmb
                i = 0
                print '>> Combining pT bins: %s' % pt_range[1]
                for ptr in pt_range[1]:
                    ptr_post = '_pt_%s' % ptr
                    for chn in channels:
                        seldir = '%s_%s_%s%s' % (region, chn, sel, ptr_post)
                        node = hists[chn][seldir][var]
                        node['data_sub'] = (node['data_obs'] - ((node['Total_R'] * bkg_sub) + node['Total_E']))

                        for path, key, val in hists[chn][seldir][var].ListObjects():
                            if i == 0:
                                hists['cmb'][seldir_cmb][var][key] = val.Clone()
                                hists['cmb'][seldir_cmb][var][key].Rebin(2)
                            else:
                                val2 = val.Clone()
                                val2.Rebin(2)
                                hists['cmb'][seldir_cmb][var][key].Add(val2)
                        i += 1

            seldir_num = '%s_sig_t%s' % (region, post)
            seldir_den = '%s_sig_l%s' % (region, post)

            true_num = '%s_iso_t_sig_t%s' % (region, post)
            true_den = '%s_iso_t_sig_l%s' % (region, post)


            fakefrac_num = hists['cmb'][seldir_num][var]['Total_J'].Clone()
            fakefrac_den = hists['cmb'][seldir_num][var]['Total_R'].Clone()
            fakefrac_den.Add(fakefrac_num)
            fakefrac_num.Divide(fakefrac_den)
            fakefrac_num.Print('range')

            results = {}
            for proc in ['data_sub', 'Total_J', 'W_H']:
                results[proc] = dict()
                res = results[proc]
                res['ratio'] = plot.MakeRatioHist(hists['cmb'][seldir_num][var][proc], hists['cmb'][seldir_den][var][proc], True, True)
                res['true_ratio'] = IntegralRatio(hists['cmb'][true_num][var][proc], hists['cmb'][true_den][var][proc])
                res['mean'] = hists['cmb'][true_den][var][proc].GetMean()
                # res['ratio'].Print('range')
                fitname = 'pol1'
                fitmin = 4
                fitmax = 20
                if region == 'endcap2':
                    fitname = 'pol1'
                    fitmin = 2
                    fitmax = 10
                res['fitres'] = res['ratio'].Fit(fitname, 'S', '', fitmin, fitmax)
                res['func'] = res['ratio'].GetListOfFunctions().At(0)
                res['err'] = ROOT.TH1D("err", "", 100, 0, 20)
                ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(res['err'], 0.68)
                res['err'].SetStats(False)
                res['err'].SetMarkerSize(0)
                res['err'].SetFillColorAlpha(2, 0.3)
                res['err'].SetLineColor(2)
                fitres = results[proc]['func'].Eval(results[proc]['mean'])
                fiterr = res['err'].GetBinError(res['err'].GetXaxis().FindFixBin(res['mean']))
                res['fit_ratio'] = (fitres, fiterr)

                plot.Set(res['ratio'], MarkerSize=0.5)
                res['ratio_cropped'] = res['ratio'].Clone()
                for ib in xrange(1, res['ratio_cropped'].GetNbinsX() + 1):
                    if res['ratio_cropped'].GetXaxis().GetBinCenter(ib) < fitmin:
                        res['ratio_cropped'].SetBinContent(ib, 0)
                        res['ratio_cropped'].SetBinError(ib, 0)

                if proc != 'data_sub':
                    plot.Set(res['ratio'], LineColor=14, MarkerColor=14, MarkerStyle=22)
                    plot.Set(res['func'], LineColor=15)
                    res['err'].SetFillColorAlpha(15, 0.2)
                    res['err'].SetLineColor(14)

            res_MC = results['Total_J']
            res_data = results['data_sub']

            # Do the fit/true ratio uncertainty
            h_ratio_num = ROOT.TH1D('h_ratio_num', '', 1, 0, 1)
            h_ratio_den = ROOT.TH1D('h_ratio_den', '', 1, 0, 1)
            h_ratio_num.SetBinContent(1, res_MC['fit_ratio'][0])
            h_ratio_num.SetBinError(1, res_MC['fit_ratio'][1])
            h_ratio_den.SetBinContent(1, res_MC['true_ratio'][0])
            h_ratio_den.SetBinError(1, res_MC['true_ratio'][1])
            h_ratio_num.Divide(h_ratio_den)
            res_MC['fit_over_true'] = [h_ratio_num.GetBinContent(1), h_ratio_num.GetBinError(1)]
            all_fit_over_true.append(res_MC['fit_over_true'])
            resline = '%-15s %-12s' % (region, pt_range[0])
            resline += ' | %.2f +/- %.2f | %.2f +/- %.2f' % (res_MC['true_ratio'][0], res_MC['true_ratio'][1], res_MC['fit_ratio'][0], res_MC['fit_ratio'][1])
            resline += ' | %.2f +/- %.2f' % (res_MC['fit_over_true'][0], res_MC['fit_over_true'][1])
            resline += ' | %.2f' % res_MC['fitres'].Prob()
            resline += ' | %.2f' % Chi2Compatibility(res_MC['fit_ratio'], res_MC['true_ratio'])[0]
            resline += ' | %.2f +/- %.2f' % (res_data['fit_ratio'][0], res_data['fit_ratio'][1])
            resline += ' | %.2f' % res_data['fitres'].Prob()

            reslines.append(resline)
            allresults[-1][-1].append(results)

            canv = ROOT.TCanvas('photon_fakes_%s_%s%s_%.1f' % (args.year, region, post, bkg_sub), 'photon_fakes_%s_%s%s_%.f' % (args.year, region, post, bkg_sub))
            pads = plot.OnePad()
            axis = results['data_sub']['ratio_cropped'].Clone()
            axis.Clear()
            axis.SetMinimum(0.0)
            axis.SetMaximum(4.0)
            plot.Set(axis.GetXaxis(), Title='Photon I_{charged} [GeV]')
            plot.Set(axis.GetYaxis(), Title='#sigma_{i#etai#eta} extrapolation')
            axis.Draw()
            results['data_sub']['ratio_cropped'].Draw('ESAME')
            results['Total_J']['ratio'].Draw('ESAME')
            results['data_sub']['err'].Draw('E3SAME')
            results['Total_J']['err'].Draw('E3SAME')

            legend = ROOT.TLegend(0.7, 0.70, 0.94, 0.93, '', 'NDC')
            legend.AddEntry(results['data_sub']['ratio_cropped'], 'Data', 'EP')
            legend.AddEntry(results['data_sub']['err'], 'Data (fit)',  'LF')
            legend.AddEntry(results['Total_J']['ratio'], 'Simulation', 'EP')
            legend.AddEntry(results['Total_J']['err'], 'Simulation (fit)', 'LF')
            legend.Draw()

            text = ROOT.TPaveText(0.50, 0.76, 0.65, 0.93, 'NDC')
            # text.AddText(args.year)
            text.AddText('p_{T} #in [%s, %s] GeV' % (pt_range[0].split('_')[0], pt_range[0].split('_')[1]))
            text.AddText('|#eta| #in [%s, %s]' % eta_ranges[region])
            text.SetTextAlign(13)
            text.SetTextFont(42)
            text.SetBorderSize(0)
            text.Draw()
            plot.DrawTitle(pads[0], args.year, 1)


            canv.Print('.pdf')
            canv.Print('.png')

    for res in reslines:
        print res

    # print len(all_fit_over_true)
    # all_fit_over_true = all_fit_over_true[:-2]
    # print len(all_fit_over_true)
    print Chi2Compatibility(*all_fit_over_true)
# The main extrap factor from fitting data
h2d = ROOT.TH2F('photon_fakes', '', 54, 30, 300, 4, array('d', [0., 1.0, 1.4442, 2.1, 2.5]))
h2d_mc_true = h2d.Clone()
h2d_mc = h2d.Clone()
h2d_stat = h2d.Clone()
h2d_bkg_syst = h2d.Clone()
h2d_const_syst = h2d.Clone()
h2d_index = h2d.Clone()

curr_index = 0

for iy in xrange(1, h2d.GetNbinsY() + 1):
    curr_index += 1
    curr_x_bin = 1
    for ix in xrange(1, h2d.GetNbinsX() + 1):
        xc = h2d.GetXaxis().GetBinCenter(ix)
        for ipt, pt_range in enumerate(config[eta_regions[iy - 1]]):
            ptmin, ptmax = [float(X) for X in pt_range[0].split('_')]
            # print ptmin, ptmax
            if xc > ptmin and xc <= ptmax:
                val = allresults[0][iy - 1][ipt]['data_sub']['fit_ratio']
                multiplier = 1.0
                if len(pt_range) > 2:
                    multiplier = pt_range[2]
                    print '>> Applying multiplier of %f' % multiplier
                val_mc_true = allresults[0][iy - 1][ipt]['Total_J']['true_ratio']
                val_mc = allresults[0][iy - 1][ipt]['Total_J']['fit_ratio']
                for h in [h2d, h2d_stat, h2d_bkg_syst, h2d_const_syst]:
                    h.SetBinContent(ix, iy, val[0] * multiplier)
                # Stat error
                err_stat = val[1] * multiplier
                err_const = common_syst * h2d.GetBinContent(ix, iy)
                err_bkg = max(abs(allresults[1][iy - 1][ipt]['data_sub']['fit_ratio'][0] - val[0]), abs(allresults[2][iy - 1][ipt]['data_sub']['fit_ratio'][0] - val[0]))
                err_tot = math.sqrt(math.pow(err_stat, 2) + math.pow(err_const, 2) + math.pow(err_bkg, 2))

                h2d_mc.SetBinContent(ix, iy, val_mc[0])
                h2d_mc.SetBinError(ix, iy, math.sqrt(math.pow(val_mc[1], 2) + math.pow(common_syst * val_mc[0], 2)))
                h2d_mc_true.SetBinContent(ix, iy, val_mc_true[0])
                h2d_mc_true.SetBinError(ix, iy, val_mc_true[1])
                h2d_stat.SetBinError(ix, iy, err_stat)
                h2d_const_syst.SetBinError(ix, iy, err_const)
                h2d_bkg_syst.SetBinError(ix, iy, err_bkg)
                h2d.SetBinError(ix, iy, err_tot)
                h2d_index.SetBinContent(ix, iy, curr_index)
                if (ix == 1) or (ipt + 1) > curr_x_bin:
                    print '%-10.0f %-10.0f %-10.2f %-10.2f %-10.2f %-10.2f' % (ptmin, ptmax, err_stat, err_const, err_bkg, err_tot)
                if (ipt + 1) > curr_x_bin:
                    curr_index += 1
                    curr_x_bin = ipt + 1
                break


fout = ROOT.TFile(args.output, 'RECREATE')


# h2d.Print('range')
# h2d_index.Print('range')

ROOT.gDirectory.WriteTObject(h2d, 'photon_fakes')
ROOT.gDirectory.WriteTObject(h2d_index, 'photon_fakes_index')
ROOT.gDirectory.WriteTObject(h2d_stat, 'photon_fakes_stat')
ROOT.gDirectory.WriteTObject(h2d_bkg_syst, 'photon_fakes_bkg_syst')
ROOT.gDirectory.WriteTObject(h2d_const_syst, 'photon_fakes_const_syst')
ROOT.gDirectory.WriteTObject(h2d_mc, 'photon_fakes_mc')
ROOT.gDirectory.WriteTObject(h2d_mc_true, 'photon_fakes_mc_true')

fout.Close()