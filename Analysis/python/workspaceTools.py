#!/usr/bin/env python
import ROOT
from array import array
import math

# Boilerplate
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.RooWorkspace.imp = getattr(ROOT.RooWorkspace, 'import')
ROOT.TH1.AddDirectory(0)
import CombineHarvester.CombineTools.plotting as plot


def TGraphAsymmErrorsToTH1D(graph):
    nbins = graph.GetN()
    bin_edges = []
    for i in xrange(0, nbins):
        bin_edges.append(graph.GetX()[i]-graph.GetEXlow()[i])
    bin_edges.append(graph.GetX()[i]+graph.GetEXhigh()[nbins-1])
    hist = ROOT.TH1D(graph.GetName(), graph.GetTitle(), nbins, array('d', bin_edges))
    for i in xrange(1, nbins+1):
        hist.SetBinContent(i, graph.GetY()[i-1])
        hist.SetBinError(i, (graph.GetEYhigh()[i-1] + graph.GetEYlow()[i-1]) / 2.)
    return hist


# Special version of the above function that infers the binning without using EXlow/high
# Instead relies on the assumption that first bin starts at zero
def TGraphAsymmErrorsToTH1DForTaus(graph):
    nbins = graph.GetN()
    bin_edges = []
    bin_edges.append(0.)
    for i in xrange(0, nbins):
        bin_edges.append(graph.GetX()[i]+(graph.GetX()[i] - bin_edges[-1]))
    print bin_edges
    hist = ROOT.TH1D(graph.GetName(), graph.GetTitle(), nbins, array('d', bin_edges))
    last_non_zero = 0.0
    for i in xrange(1, nbins+1):
        hist.SetBinContent(i, graph.GetY()[i-1])
        hist.SetBinError(i, (graph.GetEYhigh()[i-1] + graph.GetEYlow()[i-1]) / 2.)
        if graph.GetY()[i-1] != 0.:
            last_non_zero = graph.GetY()[i-1]
    # print 'Before fix'
    # hist.Print('range')
    for i in reversed(xrange(1, nbins+1)):
        if hist.GetBinContent(i) == 0.:
            hist.SetBinContent(i, last_non_zero)
        else:
            break
    # print 'After fix'
    # hist.Print('range')
    return hist

def SafeWrapHist(wsp, binvars, hist, name=None, bound=True):
    # Use the histogram name for this function unless a new name has
    # been specified
    if name is None:
        name = hist.GetName()
    # Bit of technical RooFit thing:
    # We want to use two sets of x and y variables. The first set will be
    # named specifically for the histogram in question. When we create the
    # RooDataHist RooFit will adjust the ranges and binning of these
    # RooRealVars to match the input histogram. We'll store these in
    # 'h_arglist'. The second set will contain the variables that are actually
    # used for looking up values in the RooHistFunc that will wrap around the
    # RooDataHist. These variables can be used for different RooHistFuncs with
    # different binnings and/or ranges in each, and could even be defined as
    # functions of other variables. We'll store these in 'f_arglist'.
    h_arglist = ROOT.RooArgList()
    f_arglist = ROOT.RooArgList()
    # Keep track of the relevant histogram axes (for the 1D, 2D or 3D cases)
    axes = [hist.GetXaxis()]
    if len(binvars) >= 2:
        axes.append(hist.GetYaxis())
    if len(binvars) >= 3:
        axes.append(hist.GetZaxis())
    for i, var in enumerate(binvars):
        # If the axis variable doesn't exist in the workspace already we need
        # to create it
        if not wsp.arg(var):
            # Check if the user has defined a function here
            if var.startswith('expr::'):
                funcarg = wsp.factory(var)
                var = funcarg.GetName()
            # otherwise create a normal variable
            else:
                wsp.factory('%s[0]' % var)
        # If the bound option is set we first pass the axis variables through
        # a function that adjusts the value to the xaxis minimum or maximum if
        # the user has set a value outside the axis range
        if bound:
            f_arglist.add(wsp.factory('expr::%s_range_%s("TMath::Range(%f,%f,@0)", %s)' % (
                name, var,
                axes[i].GetBinLowEdge(1),
                axes[i].GetBinUpEdge(axes[i].GetNbins()),
                var)))
        else:
            f_arglist.add(wsp.arg(var))
        # Now create the histogram-specific binning variable
        h_arglist.add(wsp.factory('%s_binningvar_%s[0]' % (name, var)))
    # First create the RooDataHist
    rdh = ROOT.RooDataHist(name, hist.GetTitle(), h_arglist, hist)
    # Then the RooHistFunc - it will create a 1->1 correspondence between the lookup variables
    # (f_arglist) and the variables that defined the histogram binning (h_arglist)
    rhf = ROOT.RooHistFunc(name, hist.GetTitle(), f_arglist, h_arglist, rdh)
    # Finally we can import it into the workspace
    wsp.imp(rhf)
    wsp.imp(hist.Clone('hist_'+name))


def MakeBinnedCategory(wsp, name, bins):
    var = wsp.factory('%s[%g,%g,%g]' % (name, bins[0], bins[0], bins[-1]))
    var.setBinning(ROOT.RooBinning(len(bins)-1, array('d', bins)))
    cat = wsp.factory('RooBinningCategory::%s_cat(%s)' % (name, name))
    return cat


def MakeBinnedCategoryFuncMap(wsp, name, bins, funcName, funcs):
    MakeBinnedCategory(wsp, name, bins)
    expr = []
    for i in xrange(len(bins)-1):
        expr.append('(@0==%i)*@%i' % (i, i+1))
    fn = wsp.factory('expr::%s("%s",%s_cat,%s)' % (funcName, '+'.join(expr), name, ','.join(funcs)))
    return fn


def ProcessDESYLeptonSFs(filename, postfix, name):
    f = ROOT.TFile(filename)
    etaBinsH = f.Get('etaBinsH')
    eta_bins = set()
    pt_bins = set()
    graphs = []
    for i in xrange(1, etaBinsH.GetNbinsX()+1):
        label = etaBinsH.GetXaxis().GetBinLabel(i)

        graph = f.Get('ZMass%s_%s' % (label, postfix))
        gr_bins = set()
        for j in xrange(0, graph.GetN()):
            gr_bins.add(graph.GetX()[j]-graph.GetEXlow()[j])
            gr_bins.add(graph.GetX()[j]+graph.GetEXhigh()[j])
        print 'Graph: %s, bins: %s' % (graph.GetName(), len(gr_bins))
        graph.Print()
        pt_bins.update(gr_bins)
        graphs.append(graph)

        eta_bins.add(etaBinsH.GetXaxis().GetBinLowEdge(i))
        eta_bins.add(etaBinsH.GetXaxis().GetBinUpEdge(i))
    result = ROOT.TH2D(name, name,
                       len(pt_bins)-1, array('d', sorted(pt_bins)),
                       len(eta_bins)-1, array('d', sorted(eta_bins)))
    for i in xrange(1, len(eta_bins)):
        for j in xrange(1, len(pt_bins)):
                result.SetBinContent(j, i, graphs[i-1].GetY()[j-1])
                result.SetBinError(j, i, (graphs[i-1].GetEYhigh()[j-1] + graphs[i-1].GetEYlow()[j-1]) / 2.)
    result.Print('range')
    return result


def HistErr(hist, err_factor=1.0, clone=True, do_ratio=True):
    if clone:
        hist_clone = hist.Clone()
    else:
        hist_clone = hist
    ndim = hist.GetDimension()
    if ndim == 1:
        for i in xrange(1, hist.GetNbinsX() + 1):
            if do_ratio:
                if hist.GetBinContent(i) > 0.:
                    hist_clone.SetBinContent(i, err_factor * hist.GetBinError(i) / hist.GetBinContent(i))
            else:
                hist_clone.SetBinContent(i, hist.GetBinContent(i) + hist.GetBinError(i) * err_factor)
    if ndim == 2:
        for i in xrange(1, hist.GetNbinsX() + 1):
            for j in xrange(1, hist.GetNbinsY() + 1):
                if do_ratio:
                    if hist.GetBinContent(i, j) > 0.:
                        hist_clone.SetBinContent(i, j, err_factor * hist.GetBinError(i, j) / hist.GetBinContent(i, j))
                else:
                    hist_clone.SetBinContent(i, j, hist.GetBinContent(i, j) + hist.GetBinError(i, j) * err_factor)
    # hist.Print('range')
    # hist_clone.Print('range')
    return hist_clone


def HistSystVariations(h_nominal, h_lo, h_hi, keepNominalErr=True, doRelative=True):
    res = h_nominal.Clone()
    for i in xrange(1, h_nominal.GetNbinsX() + 1):
        for j in xrange(1, h_nominal.GetNbinsY() + 1):
            v_nom = h_nominal.GetBinContent(i, j)
            v_err = h_nominal.GetBinError(i, j)
            v_hi = h_hi.GetBinContent(i, j)
            v_lo = h_lo.GetBinContent(i, j)
            v_av_syst = (abs(v_hi - v_nom) + abs(v_lo - v_nom)) / 2.
            if keepNominalErr:
                v_err_new = math.sqrt(math.pow(v_err, 2) + math.pow(v_av_syst, 2))
            else:
                v_err_new = v_av_syst
            res.SetBinError(i, j, v_err_new)
    return res


def ZeroErrors(h_nominal):
    res = h_nominal.Clone()
    for i in xrange(1, h_nominal.GetNbinsX() + 1):
        for j in xrange(1, h_nominal.GetNbinsY() + 1):
            res.SetBinError(i, j, 0.)
    return res




def MakeProjections(h2d_list, main_label, ix, color, marker, along='X', sublabels=['total', 'stat', 'syst']):
    res = {}
    attr = 'Projection%s' % along
    if len(h2d_list) >= 2:
        for i in xrange(len(h2d_list)):
            res[sublabels[i]] = getattr(h2d_list[i], attr)('%s_%s_proj%s_%i' % (main_label, sublabels[i], along, ix), ix, ix)

        plot.Set(res[sublabels[0]], MarkerSize=0, LineWidth=2, LineColor=color)
        plot.Set(res[sublabels[1]], LineColorAlpha=(color, 0.3), MarkerColor=color, LineWidth=8, MarkerSize=0)
        # plot.Set(res[sublabels[2]], LineColor=color, MarkerColor=color, LineWidth=2, MarkerSize=0.7, MarkerStyle=marker)
        if len(h2d_list) >= 3:
            plot.Set(res[sublabels[2]], FillColorAlpha=(0, 0), LineColor=color, MarkerColor=color, LineWidth=2, MarkerSize=0.0, MarkerStyle=marker)
    else:
        res['stat'] = getattr(h2d_list[0], attr)('%s_stat_proj%s_%i' % (main_label, along, ix), ix, ix)
        plot.Set(res['stat'], LineColor=color, MarkerColor=color, LineWidth=2, MarkerSize=0.7, MarkerStyle=marker)
    return res


def SummaryPlotsPhotonFakes(cfg):
    h_ref = cfg['h_ref'].Clone()
    main_label = cfg['main_label']
    ref_axis = h_ref.GetYaxis() if cfg['proj'] == 'X' else h_ref.GetXaxis()
    for ix in xrange(1, ref_axis.GetNbins() + 1):
        bin_label = '%s #in [%g, %g]' % (cfg['y_label'], ref_axis.GetBinLowEdge(ix), ref_axis.GetBinUpEdge(ix))
        canv = ROOT.TCanvas('%s_%s_%i' % (main_label, cfg['proj'], ix), '%s_%s_%i' % (main_label, cfg['proj'], ix))
        pads = plot.OnePad()
        pads[0].cd()

        text = ROOT.TPaveText(0.17, 0.84, 0.6, 0.93, 'NDC')
        legend = ROOT.TLegend(0.6, 0.75, 0.94, 0.93, '', 'NDC')
        data_hists = MakeProjections(cfg['data'], '%s_data' % main_label, ix, color=cfg.get('data_colour', 4), marker=21, along=cfg['proj'], sublabels=cfg['data_labels'])
        # mc_hists = MakeProjections(cfg['mc'], '%s_mc' % main_label, ix, color=cfg.get('mc_colour', 2), marker=20, along=cfg['proj'])
        # ratio_hists = MakeProjections(cfg['ratio'], '%s_ratio' % main_label, ix, color=cfg.get('ratio_colour', 1), marker=21, along=cfg['proj'])

        if 'const_syst' in data_hists:
            data_hists['total_err'] = data_hists['total'].Clone()
            plot.Set(data_hists['total_err'], FillColorAlpha=(data_hists['total_err'].GetLineColor(), 0.3), MarkerColor=data_hists['total_err'].GetLineColor(), LineWidth=0, MarkerSize=0)
            plot.Set(data_hists['bkg_syst'], FillColorAlpha=(0, 0), LineColor=2, LineWidth=2, MarkerSize=0)
            plot.Set(data_hists['const_syst'], FillColorAlpha=(0, 0), LineColor=1, LineWidth=2, MarkerSize=0)
            plot.Set(data_hists['stat'], FillColorAlpha=(0, 0), LineColor=4, LineWidth=2, LineStyle=2)
            for h in ['stat', 'const_syst', 'bkg_syst']:
                data_hists[h + '_hi'] = data_hists[h].Clone()
                data_hists[h + '_lo'] = data_hists[h].Clone()
            for ib in xrange(1, data_hists['const_syst'].GetNbinsX() + 1):
                for h in ['stat', 'const_syst', 'bkg_syst']:
                    data_hists[h + '_hi'].SetBinContent(ib, data_hists[h].GetBinContent(ib) + data_hists[h].GetBinError(ib))
                    data_hists[h + '_lo'].SetBinContent(ib, data_hists[h].GetBinContent(ib) - data_hists[h].GetBinError(ib))
                    data_hists[h + '_hi'].SetBinError(ib, 1E-6)
                    data_hists[h + '_lo'].SetBinError(ib, 1E-6)
                data_hists['total'].SetBinError(ib, 1E-6)
                data_hists['total'].SetBinError(ib, 1E-6)

            data_hists['total'].Draw('ESAME')
            data_hists['total_err'].Draw('E2SAME')
            data_hists['const_syst_hi'].Draw('ESAME')
            data_hists['const_syst_lo'].Draw('ESAME')
            data_hists['bkg_syst_hi'].Draw('ESAME')
            data_hists['bkg_syst_lo'].Draw('ESAME')
            data_hists['stat_hi'].Draw('ESAME')
            data_hists['stat_lo'].Draw('ESAME')

            legend.AddEntry(data_hists['total'], 'Data', 'L')
            legend.AddEntry(data_hists['stat'], '  Statistical', 'L')
            legend.AddEntry(data_hists['const_syst'], '  Non-closure syst.', 'L')
            legend.AddEntry(data_hists['bkg_syst'], '  Prompt subtraction syst.', 'L')
            legend.AddEntry(data_hists['total_err'], '  Total', 'F')
            # legend.AddEntry(mc_hists['stat'], 'Simulation', 'P')
            # legend.AddEntry(mc_hists['stat'], '  Statistical', 'E')
        elif 'mc' in data_hists:
            plot.Set(data_hists['mc'], LineColor=13, LineWidth=2, MarkerSize=0)
            plot.Set(data_hists['mc_true'], LineColor=9, LineWidth=2, MarkerSize=0)
            data_hists['mc_err'] = data_hists['mc'].Clone()
            plot.Set(data_hists['mc_err'], FillColorAlpha=(data_hists['mc_err'].GetLineColor(), 0.3), LineWidth=0, MarkerSize=0)
            data_hists['mc_true_err'] = data_hists['mc_true'].Clone()
            plot.Set(data_hists['mc_true_err'], FillColorAlpha=(data_hists['mc_true_err'].GetLineColor(), 0.3), LineWidth=0, MarkerSize=0)
            for ib in xrange(1, data_hists['mc'].GetNbinsX() + 1):
                data_hists['mc'].SetBinError(ib, 1E-6)
                data_hists['mc_true'].SetBinError(ib, 1E-6)

            data_hists['mc'].Draw('ESAME')
            data_hists['mc_err'].Draw('E2SAME')
            data_hists['mc_true'].Draw('ESAME')
            data_hists['mc_true_err'].Draw('E2SAME')

            legend.AddEntry(data_hists['mc'], 'Simulation - fit', 'LF')
            legend.AddEntry(data_hists['mc_true'], 'Simulation - truth', 'LF')
            # legend.AddEntry(mc_hists['stat'], 'Simulation', 'PE')

        axis = plot.GetAxisHist(pads[0])
        plot.Set(axis, Minimum=cfg['y_range'][0], Maximum=cfg['y_range'][1])

        plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.40)

        axis.GetYaxis().SetTitle('#sigma_{i#etai#eta} extrapolation')
        axis.GetXaxis().SetTitle(cfg['x_axis_title'])


        pads[0].cd()
        legend.Draw()
        text.AddText(cfg['main_text'])
        text.AddText(bin_label)
        text.SetTextAlign(13)
        text.SetBorderSize(0)
        text.Draw()
        pads[0].SetGrid(1, 1)
        if cfg['logx']:
            pads[0].SetLogx(True)
        canv.Print('.png')
        canv.Print('.pdf')


def SummaryPlots(cfg):
    h_ref = cfg['h_ref'].Clone()
    main_label = cfg['main_label']
    ref_axis = h_ref.GetYaxis() if cfg['proj'] == 'X' else h_ref.GetXaxis()
    for ix in xrange(1, ref_axis.GetNbins() + 1):
        bin_label = '%s #in [%g, %g]' % (cfg['y_label'], ref_axis.GetBinLowEdge(ix), ref_axis.GetBinUpEdge(ix))
        canv = ROOT.TCanvas('%s_%s_%i' % (main_label, cfg['proj'], ix), '%s_%s_%i' % (main_label, cfg['proj'], ix))
        pads = plot.TwoPadSplit(cfg.get('ratio_split', 0.4), 0.01, 0.01)
        pads[0].cd()

        text = ROOT.TPaveText(0.17, 0.84, 0.6, 0.93, 'NDC')
        legend = ROOT.TLegend(0.6, 0.75, 0.94, 0.93, '', 'NDC')
        data_hists = MakeProjections(cfg['data'], '%s_data' % main_label, ix, color=cfg.get('data_colour', 4), marker=21, along=cfg['proj'])
        mc_hists = MakeProjections(cfg['mc'], '%s_mc' % main_label, ix, color=cfg.get('mc_colour', 2), marker=20, along=cfg['proj'])
        ratio_hists = MakeProjections(cfg['ratio'], '%s_ratio' % main_label, ix, color=cfg.get('ratio_colour', 1), marker=21, along=cfg['proj'])

        if 'total' in data_hists:
            data_hists['total'].Draw('E2SAME')
        if 'syst' in data_hists:
            data_hists['syst'].Draw('E0X0SAME')
        if 'stat' in data_hists:
            data_hists['stat'].Draw('E1X0SAME')
        if 'stat' in mc_hists:
            mc_hists['stat'].Draw('E1X0PSAME')

        if 'syst' in data_hists:
            legend.AddEntry(data_hists['stat'], 'Data', 'P')
            legend.AddEntry(data_hists['stat'], '  Statistical', 'E')
            legend.AddEntry(data_hists['syst'], '  Systematic', 'E')
            legend.AddEntry(data_hists['total'], '  Total', 'F')
            legend.AddEntry(mc_hists['stat'], 'Simulation', 'P')
            legend.AddEntry(mc_hists['stat'], '  Statistical', 'E')
        else:
            legend.AddEntry(data_hists['stat'], 'Data', 'PE')
            legend.AddEntry(mc_hists['stat'], 'Simulation', 'PE')

        axis = plot.GetAxisHist(pads[0])
        plot.Set(axis, Minimum=cfg['y_range'][0], Maximum=cfg['y_range'][1])

        plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.40)

        axis.GetYaxis().SetTitle('Efficiency')

        pads[1].cd()

        if 'total' in ratio_hists:
            ratio_hists['total'].Draw('E2SAME')
        if 'syst' in ratio_hists:
            ratio_hists['syst'].Draw('E0X0SAME')
        if 'stat' in ratio_hists:
            ratio_hists['stat'].Draw('E1X0SAME')
        plot.SetupTwoPadSplitAsRatio(pads,
                                     plot.GetAxisHist(pads[0]),
                                     plot.GetAxisHist(pads[1]),
                                     'Data/Sim', True, cfg['ratio_range'][0], cfg['ratio_range'][1])
        r_axis = plot.GetAxisHist(pads[1])

        if cfg['logx']:
            pads[1].SetLogx(True)
            r_axis.GetXaxis().SetMoreLogLabels(True)
            r_axis.GetXaxis().SetNoExponent(True)
        r_axis.GetXaxis().SetTitle(cfg['x_axis_title'])
        r_axis.GetXaxis().SetTitleOffset(ROOT.gStyle.GetTitleXOffset())
        pads[1].SetGrid(1, 1)
        pads[1].RedrawAxis('g')

        pads[0].cd()
        legend.Draw()
        text.AddText(cfg['main_text'])
        text.AddText(bin_label)
        text.SetTextAlign(13)
        text.SetBorderSize(0)
        text.Draw()
        pads[0].SetGrid(1, 1)
        if cfg['logx']:
            pads[0].SetLogx(True)
        canv.Print('.png')
        canv.Print('.pdf')