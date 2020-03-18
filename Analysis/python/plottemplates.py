import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
import ROOT
from copy import deepcopy
from array import array

ROOT.TH1.SetDefaultSumw2(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()

AUTO_COLOURS = [
    [95, 173, 86],
    [242, 193, 78],
    [247, 129, 84],
    [81, 163, 163],
    [117, 72, 94],
    [19, 41, 61],
    [237, 37, 78]
]

def HistSum(hdict, label_list):
    return sum([hdict[X] for X in label_list[1:]], hdict[label_list[0]])


def VariableRebin(hist, binning):
    newhist = hist.Rebin(len(binning) - 1, "", array('d', binning))
    return newhist


def UpdateDict(current, update):
    cfg = deepcopy(current)
    cfg.update(update)
    return cfg


DEFAULT_CFG = {
    # full filename will be [outdir]/[prefix]name[postfix].[ext]:
    'type': 'datamc',           # or multihist
    'pads': None,               # externally sourced TPads
    'outdir': '',               # output directory
    'prefix': '',               # filename prefix
    'postfix': '',              # filename postfix
    'logx': False,              # Draw x-axis in log-scale
    'logy': False,              # Draw y-axis in log-scale
    'logy_min': 1E-3,
    'ratio': True,             # Draw the ratio plot?
    'fraction': False,             # Draw the ratio plot?
    'purity': False,
    'ratio_pad_frac': 0.27,
    'ratio_y_range': [0.61, 1.39],  # Range of the ratio y-axis
    'ratio_y_title': 'Obs/Exp',
    'x_range': [],                  # Restrict the x-axis range shown
    'rebin': 0,                     # Rebin by this factor
    'rebinvar': [],                 # Rebin to this list of bin edges
    'x_title': '',               # Title on the x-axis
    'y_title': 'Events',         # Title on the y-axis
    'divwidth': True,           # Divide all histogram contents by their bin widths
    'layout': 'data_fakes',
    'legend_pos': [0.67, 0.66, 0.90, 0.91], # Legend position in NDC (x1, y1, x2, y2)
    'legend_padding': 0.05,      # Automatically increase the y-axis range to ensure the legend does not overlap the histograms (argument is fraction of frame height to pad)
    'legend_show_yields': False,  # Add the yields to the legends
    'data_name': 'data_obs',     # Name of the TH1 to take for data
    'main_logo': 'CMS',
    'sub_logo': 'Internal',
    'top_title_right': '35.9 fb^{-1} (13 TeV)',
    'top_title_left': '',
    'hide_data': False,
    'auto_top_title_right': True,
    'global_hist_opts': {
        'draw_opts': 'E0',
        'legend_opts': 'L',
        'marker_size': 0.6,
        'line_width': 3
    },
    'norm_to': 0,
    'overlays': [
        # {
        #     "name": "systUp",
        #     "entries": "total",
        #     "hist_postfix": "_CMS_scale_met_jesUp",
        #     "legend": "systUp",
        #     "color": 2
        # },
        # {
        #     "name": "systDown",
        #     "entries": "total",
        #     "hist_postfix": "_CMS_scale_met_jesDown",
        #     "legend": "systDown",
        #     "color": 4
        # }
    ]
}


def MakeMultiHistPlot(name, outdir, hists, cfg, layout, ratios=None):
    copyhists = {}
    for hname, h in hists.iteritems():
        if len(cfg['rebinvar']):
            copyhists[hname] = VariableRebin(h, cfg['rebinvar'])
        else:
            copyhists[hname] = h.Clone()
        if cfg['norm_to'] > 0.:
            copyhists[hname].Scale(cfg['norm_to'] / copyhists[hname].Integral())
        if cfg['divwidth']:
            copyhists[hname].Scale(1., 'width')

    hists = copyhists

    # Canvas and pads
    if cfg['pads'] is not None:
        pads = cfg['pads']
    elif cfg['ratio'] or cfg['fraction']:
        canv = ROOT.TCanvas(name, name)
        if cfg['purity']:
            pads = plot.MultiRatioSplit([0.27, 0.13], [0.005, 0.005], [0.005, 0.005])
        else:
            pads = plot.TwoPadSplit(cfg['ratio_pad_frac'], 0.01, 0.01)
    else:
        canv = ROOT.TCanvas(name, name)
        pads = plot.OnePad()

    # Allow the user to skip specifying a list of entries for a given plot element.
    # If it's not been specified then we will add it manually
    for info in layout:
        if 'entries' not in info:
            info['entries'] = [info['name']]
        for opt in cfg['global_hist_opts']:
            if opt not in info:
                info[opt] = cfg['global_hist_opts'][opt]

    h_data = None
    if cfg['type'] == 'datamc':
        # Get the data and create axis hist
        h_data = hists[cfg['data_name']]
    else:
        h_data = hists[layout[0]['entries'][0]]

    if isinstance(h_data, ROOT.TH2):
        print 'TH2: aborting!'
        return

    h_axes = [h_data.Clone() for x in pads]
    for h in h_axes:
        if len(cfg['x_range']):
            h.GetXaxis().SetRangeUser(*cfg['x_range'])
        h.Reset()

    build_h_tot = True
    h_tot = None
    if 'TotalProcs' in hists:
        h_tot = hists['TotalProcs']
        build_h_tot = False
    if cfg['type'] != 'datamc':
        build_h_tot = False

    if isinstance(cfg['x_title'], list) or isinstance(cfg['x_title'], tuple):
        x_title = cfg['x_title'][0]
        units = cfg['x_title'][1]
    else:
        x_title = cfg['x_title']
        units = ''

    if x_title == '' and h_data.GetXaxis().GetTitle() != '':
        x_title = h_data.GetXaxis().GetTitle()

    if ':' in x_title:
        units = x_title.split(':')[1]
        x_title = x_title.split(':')[0]

    if cfg['logy']:
        pads[0].SetLogy()
        h_axes[0].SetMinimum(cfg['logy_min'])

    rpad_idx = len(pads) - 1

    if cfg['ratio'] or cfg['fraction']:
        plot.StandardAxes(h_axes[rpad_idx].GetXaxis(), h_axes[0].GetYaxis(), x_title, units)
    else:
        plot.StandardAxes(h_axes[0].GetXaxis(), h_axes[0].GetYaxis(), x_title, units)
    h_axes[0].Draw()

    # A dict to keep track of the hists
    h_store = {}
    p_store = {}

    legend = ROOT.TLegend(*(cfg['legend_pos'] + ['', 'NBNDC']))
    stack = ROOT.THStack()
    purity_stack = ROOT.THStack()

    curr_auto_colour = 0
    for info in layout:
        hist = hists[info['entries'][0]]
        if 'color' in info:
            col = info['color']
        else:
            col = AUTO_COLOURS[curr_auto_colour]
            if curr_auto_colour == (len(AUTO_COLOURS) - 1):
                curr_auto_colour = 0
            else:
                curr_auto_colour += 1
        # col = info['color']
        if isinstance(col, list):
            col = ROOT.TColor.GetColor(*col)
        # print info['line_width']
        if cfg['type'] == 'multihist':
            plot.Set(hist, FillColor=col, MarkerColor=col, LineColor=col, Title=info['legend'], MarkerSize=info['marker_size'], LineWidth=info['line_width'])
        else:
            plot.Set(hist, FillColor=col, MarkerColor=col, Title=info['legend'], MarkerSize=info['marker_size'], LineWidth=info['line_width'])
        if len(info['entries']) > 1:
            for other in info['entries'][1:]:
                hist.Add(hists[other])
        h_store[info['name']] = hist
        p_store[info['name']] = hist.Clone()
        if build_h_tot:
            if h_tot is None:
                h_tot = hist.Clone()
            else:
                h_tot.Add(hist)
        if cfg['type'] == 'datamc':
            stack.Add(hist)
        else:
            hist.Draw('SAME%s' % info['draw_opts'])

    # h_tot_purity = h_tot.Clone()
    for info in layout:
        p_store[info['name']].Divide(h_tot)
        purity_stack.Add(p_store[info['name']])

    if cfg['type'] == 'datamc':
        h_tot.SetFillColor(plot.CreateTransparentColor(12, 0.3))
        h_tot.SetMarkerSize(0)
        legend.AddEntry(h_data, 'Observed', 'PL')

    # Build overlays
    for info in cfg['overlays']:
        hist = None
        input_list = []
        if isinstance(info['entries'], str):
            input_list = list(all_input_hists)
        else:
            input_list = list(info['entries'])
        updated_list = []
        for xh in input_list:
            if xh + info['hist_postfix'] in hists:
                updated_list.append(xh + info['hist_postfix'])
            else:
                updated_list.append(xh)
        print updated_list
        hist = HistSum(hists, updated_list)
        col = info['color']
        if isinstance(col, list):
            col = ROOT.TColor.GetColor(*col)
        plot.Set(hist, LineColor=col, LineWidth=1, MarkerSize=0, Title=info['legend'])
        for ib in xrange(1, hist.GetNbinsX() + 1):
            hist.SetBinError(ib, 1E-7)
        h_store[info['name']] = hist

    if cfg['type'] == 'datamc':
        for ele in reversed(layout):
            legend.AddEntry(h_store[ele['name']], ele['legend'], ele['legend_opts'])
    else:
        for ele in layout:
            leg_extra = ''
            if cfg['legend_show_yields']:
                leg_extra = ' (%.1f)' % h_store[ele['name']].Integral('width' if cfg['divwidth'] else '')
            legend.AddEntry(h_store[ele['name']], ele['legend'] + leg_extra, ele['legend_opts'])

    if cfg['type'] == 'datamc':
        bkg_uncert_label = 'Stat. Uncertainty'
        if not build_h_tot:
            bkg_uncert_label = 'Uncertainty'
        legend.AddEntry(h_tot, bkg_uncert_label, 'F')

        stack.Draw('HISTSAME')
        h_tot.Draw("E2SAME")

        for info in cfg['overlays']:
            h_store[info['name']].Draw('HISTSAME')

        if not cfg['hide_data']:
            h_data.Draw('E0SAME')

    for info in cfg['overlays']:
        legend.AddEntry(h_store[info['name']], info['legend'], 'L')


    plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.35)
    legend.Draw()
    # if cfg['legend_padding'] > 0.:
    #     plot.FixBoxPadding(pads[0], legend, cfg['legend_padding'])

    # Do the ratio plot
    r_store = {}
    r_data = None
    r_tot = None
    if cfg['ratio'] or cfg['fraction']:
        pads[rpad_idx].cd()
        pads[rpad_idx].SetGrid(0, 1)
        h_axes[rpad_idx].Draw()

        if cfg['type'] == 'datamc' and cfg['ratio']:
            r_data = plot.MakeRatioHist(h_data, h_tot, True, False)
            r_tot = plot.MakeRatioHist(h_tot, h_tot, True, False)
            r_tot.Draw('E2SAME')
            for info in cfg['overlays']:
                r_store[info['name']] = plot.MakeRatioHist(h_store[info['name']], h_tot, False, False)
                r_store[info['name']].Draw('SAME')
            if not cfg['hide_data']:
                r_data.Draw('SAME')

        if cfg['type'] == 'datamc' and cfg['fraction']:
            r_frac = plot.MakeRatioHist(h_tot, h_data, True, True)
            r_frac.Draw('SAME')
            plot.SetupTwoPadSplitAsRatio(
                pads, plot.GetAxisHist(
                    pads[0]), plot.GetAxisHist(pads[rpad_idx]), 'Exp/Obs', True, 0.0, 0.5)

        if ratios is not None:
            for info in ratios:
                if 'type' in info and info['type'] == 'binomial':
                    rhist = h_store[info['num']].Clone()
                    rhist.Divide(h_store[info['num']], h_store[info['den']], 1., 1., "B")
                elif 'type' in info and info['type'] == 'noerror':
                    rhist = plot.MakeRatioHist(h_store[info['num']], h_store[info['den']], True, False)
                else:
                    rhist = plot.MakeRatioHist(h_store[info['num']], h_store[info['den']], True, True)
                r_store['%s_%s' % (info['num'], info['den'])] = rhist
                rhist.Draw('SAMEE0')

        plot.SetupTwoPadSplitAsRatio(
            pads, plot.GetAxisHist(
                pads[0]), plot.GetAxisHist(pads[rpad_idx]), cfg['ratio_y_title'], True, *(cfg['ratio_y_range']))

    if cfg['purity']:
        pads[1].cd()
        h_axes[1].Draw()
        plot.SetupTwoPadSplitAsRatio(pads, plot.GetAxisHist(
                    pads[0]), plot.GetAxisHist(pads[1]), 'f', True, 0, 1)
        plot.Set(h_axes[1].GetXaxis(), TitleSize=0, LabelSize=0)
        plot.Set(h_axes[1].GetYaxis(), Ndivisions=(502, False))
        purity_stack.Draw('HISTSAME')
        # purity_stack.Print()
        h_axes[1].SetMinimum(0)
        h_axes[1].SetMaximum(1)
        pads[1].RedrawAxis()

    # Go back and tidy up the axes and frame
    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    # CMS logo
    plot.DrawCMSLogo(pads[0], cfg['main_logo'], cfg['sub_logo'], 11, 0.045, 0.05, 1.0, '', 1.0)
    plot.DrawTitle(pads[0], cfg['top_title_left'], 1)
    if cfg['auto_top_title_right']:
        title_right = h_data.GetTitle()
        if title_right.startswith('lumi:'):
            plot.DrawTitle(pads[0], title_right.replace('lumi:', ''), 3)

    latex = ROOT.TLatex()
    plot.Set(latex, NDC=None, TextFont=42, TextSize=0.03)
    # latex.DrawLatex(0.20, 0.75, args.title)
    # plot.DrawTitle(pads[0], args.title, 1)

    # ... and we're done
    if cfg['pads'] is None:
        canv.Print(outdir + '/' + cfg['prefix'] + name + cfg['postfix'] + '.png')
        canv.Print(outdir + '/' + cfg['prefix'] + name + cfg['postfix'] + '.pdf')

    outobjs = {}
    outobjs['axes'] = h_axes
    outobjs['hists'] = hists
    outobjs['stack'] = stack
    outobjs['purity_stack'] = purity_stack
    outobjs['h_tot'] = h_tot
    outobjs['legend'] = legend
    outobjs['r_store'] = r_store
    outobjs['p_store'] = p_store
    outobjs['r_data'] = r_data
    outobjs['r_tot'] = r_tot
    return outobjs

# def MakeDataMCPlot(name, outdir, hists, cfg, layouts):
#     copyhists = {}
#     for hname, h in hists.iteritems():
#         if len(cfg['rebinvar']):
#             copyhists[hname] = VariableRebin(h, cfg['rebinvar'])
#         else:
#             copyhists[hname] = h.Clone()
#         if cfg['divwidth']:
#             copyhists[hname].Scale(1., 'width')

#     hists = copyhists

#     # Canvas and pads
#     canv = ROOT.TCanvas(name, name)
#     if cfg['ratio']:
#         pads = plot.TwoPadSplit(0.27, 0.01, 0.01)
#     else:
#         pads = plot.OnePad()

#     # Get the data and create axis hist
#     h_data = hists[cfg['data_name']]
#     # h_data = Getter(file, '%s/data_obs' % target, True)
#     if isinstance(h_data, ROOT.TH2):
#         print 'TH2: aborting!'
#         return

#     h_axes = [h_data.Clone() for x in pads]
#     for h in h_axes:
#         if len(cfg['x_range']):
#             h.GetXaxis().SetRangeUser(*cfg['x_range'])
#         h.Reset()

#     build_h_tot = True
#     h_tot = None
#     if 'TotalProcs' in hists:
#         h_tot = hists['TotalProcs']
#         build_h_tot = False

#     x_title = cfg['x_title'][0]
#     units = cfg['x_title'][1]

#     if x_title == '' and h_data.GetXaxis().GetTitle() != '':
#         x_title = h_data.GetXaxis().GetTitle()

#     if ':' in x_title:
#         units = x_title.split(':')[1]
#         x_title = x_title.split(':')[0]

#     if cfg['logy']:
#         pads[0].SetLogy()
#         h_axes[0].SetMinimum(0.001)

#     if cfg['ratio']:
#         plot.StandardAxes(h_axes[1].GetXaxis(), h_axes[0].GetYaxis(), x_title, units)
#     else:
#         plot.StandardAxes(h_axes[0].GetXaxis(), h_axes[0].GetYaxis(), x_title, units)
#     h_axes[0].Draw()

#     # A dict to keep track of the hists
#     h_store = {}

#     layout = layouts[cfg['layout']]

#     stack = ROOT.THStack()
#     legend = ROOT.TLegend(*(cfg['legend_pos'] + ['', 'NBNDC']))

#     for info in layout:
#         hist = hists[info['entries'][0]]
#         col = info['color']
#         if isinstance(col, list):
#             col = ROOT.TColor.GetColor(*col)
#         plot.Set(hist, FillColor=col, Title=info['legend'])
#         if len(info['entries']) > 1:
#             for other in info['entries'][1:]:
#                 hist.Add(hists[other])
#         h_store[info['name']] = hist
#         if build_h_tot:
#             if h_tot is None:
#                 h_tot = hist.Clone()
#             else:
#                 h_tot.Add(hist)
#         stack.Add(hist)

#     h_tot.SetFillColor(plot.CreateTransparentColor(12, 0.3))
#     h_tot.SetMarkerSize(0)

#     legend.AddEntry(h_data, 'Observed', 'PL')
#     for ele in reversed(layout):
#         legend.AddEntry(h_store[ele['name']], '', 'F')
#     bkg_uncert_label = 'Stat. Uncertainty'
#     if not build_h_tot:
#         bkg_uncert_label = 'Uncertainty'
#     legend.AddEntry(h_tot, bkg_uncert_label, 'F')

#     stack.Draw('HISTSAME')
#     h_tot.Draw("E2SAME")
#     if not cfg['hide_data']:
#         h_data.Draw('E0SAME')

#     plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.35)
#     legend.Draw()
#     if cfg['legend_padding'] > 0.:
#         plot.FixBoxPadding(pads[0], legend, cfg['legend_padding'])

#     # Do the ratio plot
#     if cfg['ratio']:
#         pads[1].cd()
#         pads[1].SetGrid(0, 1)
#         h_axes[1].Draw()

#         r_data = plot.MakeRatioHist(h_data, h_tot, True, False)
#         r_tot = plot.MakeRatioHist(h_tot, h_tot, True, False)
#         r_tot.Draw('E2SAME')
#         if not cfg['hide_data']:
#             r_data.Draw('SAME')

#         plot.SetupTwoPadSplitAsRatio(
#             pads, plot.GetAxisHist(
#                 pads[0]), plot.GetAxisHist(pads[1]), 'Obs/Exp', True, 0.61, 1.69)

#     # Go back and tidy up the axes and frame
#     pads[0].cd()
#     pads[0].GetFrame().Draw()
#     pads[0].RedrawAxis()

#     # CMS logo
#     plot.DrawCMSLogo(pads[0], cfg['main_logo'], cfg['sub_logo'], 11, 0.045, 0.05, 1.0, '', 1.0)
#     plot.DrawTitle(pads[0], cfg['top_title_left'], 1)
#     plot.DrawTitle(pads[0], cfg['top_title_right'], 3)

#     latex = ROOT.TLatex()
#     plot.Set(latex, NDC=None, TextFont=42, TextSize=0.03)
#     latex.DrawLatex(0.20, 0.75, args.title)
#     # plot.DrawTitle(pads[0], args.title, 1)

#     # ... and we're done
#     canv.Print(outdir + '/' + cfg['prefix'] + name + cfg['postfix'] + '.png')
#     canv.Print(outdir + '/' + cfg['prefix'] + name + cfg['postfix'] + '.pdf')

