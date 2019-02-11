import ROOT
import argparse
import json
import itertools
from copy import deepcopy
from pprint import pprint
import CombineHarvester.CombineTools.plotting as plot
from Acorn.Analysis.analysis import *
from array import array
from Acorn.Analysis.plottemplates import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle(width=900, height=900)


def NormTH2InColumns(h):
    for i in xrange(1, h.GetNbinsX() + 1):
        col_tot = 0.
        for j in xrange(1, h.GetNbinsY() + 1):
            col_tot += h.GetBinContent(i, j)
        if col_tot == 0.:
            continue
        for j in xrange(1, h.GetNbinsY() + 1):
            h.SetBinContent(i, j, h.GetBinContent(i, j) / col_tot)


parser = argparse.ArgumentParser()
parser.add_argument('--load', action='store_true')
parser.add_argument('--years', default='2016,2017,2018')
parser.add_argument('--input', default='root://eoscms.cern.ch//store/cmst3/user/agilbert/190208')
parser.add_argument('--output', default='output/plots')
args = parser.parse_args()

years = args.years.split(',')
inprefix = args.input
outdir = args.output
# inprefix = '/home/files'
indirs = {
    '2016': inprefix + '/wgamma_2016_v3/WGamma_',
    '2017': inprefix + '/wgamma_2017_v3/WGamma_',
    '2018': inprefix + '/wgamma_2018_v3/WGamma_',
}

samples = {
    '2016': {
        'data_obs': 'SingleMuon',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG-inc': 'WGToLNuG-madgraphMLM',
        'WG-NLO': 'WGToLNuG-amcatnloFXFX'
    },
    '2017': {
        'data_obs': 'SingleMuon',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG-inc': 'WGToLNuG-madgraphMLM',
        'WG-NLO': 'WGToLNuG-madgraphMLM-stitched'  # No NLO sample yet
    },
    '2018': {
        'data_obs': 'SingleMuon',
        'WG': 'WGToLNuG-madgraphMLM-stitched',
        'WG-inc': 'WGToLNuG-madgraphMLM',
        'WG-NLO': 'WGToLNuG-amcatnloFXFX'
    }
}

sample_files = {}
sample_cfgs = {}
lumi_fb = {}

for year in years:
    sample_files[year] = {
        key: indirs[year] + value + '.root' for (key, value) in samples[year].iteritems()
    }

    with open('input/cfg_wgamma_%s_v3.json' % year) as jsonfile:
        cfg = json.load(jsonfile)
        sample_cfgs[year] = cfg['samples']

    for sample in samples[year]:
        if sample == 'data_obs':
            lumi_fb[year] = sample_cfgs[year][samples[year][sample]]['lumi'] / 1000.
            continue
        f = ROOT.TFile.Open(sample_files[year][sample])
        sample_cfgs[year][samples[year][sample]]['events'] = f.Get('counters').GetBinContent(2)
        f.Close()

hists = Node()

if not args.load:
    list_of_vars = [
        ('lhe_p0_pt_wide', 'lhe_p0_pt', [0, 10, 20, 30, 40, 50, 100, 150, 210, 300, 420, 600, 850, 1200]),
        ('lhe_p0_pt', 'lhe_p0_pt', (50, 0, 1000)),
        ('p0_pt', 'gen_p0_pt', [0, 10, 20, 30, 40, 50, 100, 150, 210, 300, 420, 600, 850, 1200]),
        ('p0_eta', 'gen_p0_eta', (20, -3, 3)),
        ('l0_pt', 'gen_l0_pt', [0, 10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500]),
        ('l0_eta', 'gen_l0_eta', (20, -3, 3)),
    ]

    for label, var, binning in list_of_vars:
        X = SelectionManager()
        X['baseline_m'] = 'gen_pdgid==13 && is_wg_gen && gen_l0_pt > 30 && abs(gen_l0_eta) < 2.4 && gen_p0_pt > 30 && abs(gen_p0_eta) < 2.5 && (abs(gen_p0_eta) < 1.4442 || abs(gen_p0_eta) > 1.566) && gen_l0p0_dr > 0.7'
        X['baseline_e'] = 'gen_pdgid==11 && is_wg_gen && gen_l0_pt > 35 && abs(gen_l0_eta) < 2.4 && (abs(gen_l0_eta) < 1.4442 || abs(gen_l0_eta) > 1.566) && gen_p0_pt > 30 && abs(gen_p0_eta) < 2.5 && (abs(gen_p0_eta) < 1.4442 || abs(gen_p0_eta) > 1.566) && gen_l0p0_dr > 0.7'
        X['muon_id'] = '$baseline_m && n_pre_m>=1 && l0_pdgid == 13'
        X['elec_id'] = '$baseline_e && n_pre_e>=1 && l0_pdgid == 11'
        X['elec_match'] = '$elec_id && gen_l0_match'
        X['muon_iso'] = '$muon_id && l0_iso < 0.15'
        X['elec_iso'] = '$elec_id'
        X['muon_trg'] = '$muon_iso && (l0_trg || (l0_pt > 55 && l0_trg_2))'
        X['elec_trg'] = '$elec_iso && l0_trg'
        X['muon_trg_2'] = '$muon_iso && l0_trg'
        X['elec_trg_2'] = '$elec_iso && l0_trg'
        X['photon_id_m'] = '$muon_trg_2 && n_pre_p>=1 && p0_medium'
        X['photon_pix_m'] = '$photon_id_m && !p0_haspix'
        X['photon_id_e'] = '$elec_trg_2 && n_pre_p>=1 && p0_medium'
        X['photon_match_e'] = '$photon_id_e && gen_p0_match'
        X['photon_pix_e'] = '$photon_match_e && !p0_haspix'
        X['eft_m'] = '$photon_pix_m && l0_pt>80 && gen_met>80 && p0_pt>150 && l0p0_dr>3.0'
        X['eft_e'] = '$photon_pix_e && l0_pt>80 && gen_met>80 && p0_pt>150 && l0p0_dr>3.0'

# 'm0_pt>80 && met>80 && p0_pt>150 && m0p0_dr>3.0'
        for year in years:
            for sel in X.storage.keys():
                for sa in ['WG', 'WG-inc', 'WG-NLO']:
                    hists[year][label][sa][sel] = Hist('TH1D', sample=sa, var=[var], binning=binning, sel=X.get('$' + sel), wt='wt_pu*wt_def')
            for sel in ['photon_pix_m', 'photon_pix_e', 'eft_m', 'eft_e']:
                for sa in ['WG']:
                    hists[year]['sphi_response'][sa][sel] = Hist('TH2F', sample=sa, var=['gen_sphi', 'reco_sphi'], binning=(12, -3.142, 3.142, 12, -3.142, 3.142), sel=X.get('$' + sel), wt='wt_pu*wt_def')
                    hists[year]['sphi_xy_response'][sa][sel] = Hist('TH2F', sample=sa, var=['gen_sphi', 'reco_xy_sphi'], binning=(12, -3.142, 3.142, 12, -3.142, 3.142), sel=X.get('$' + sel), wt='wt_pu*wt_def')
                    hists[year]['sphi_puppi_response'][sa][sel] = Hist('TH2F', sample=sa, var=['gen_sphi', 'reco_puppi_sphi'], binning=(12, -3.142, 3.142, 12, -3.142, 3.142), sel=X.get('$' + sel), wt='wt_pu*wt_def')
                    hists[year]['phi_response'][sa][sel] = Hist('TH2F', sample=sa, var=['gen_phi', 'reco_phi'], binning=(12, -3.142, 3.142, 12, -3.142, 3.142), sel=X.get('$' + sel), wt='wt_pu*wt_def')
                    hists[year]['true_gen_phi'][sa][sel] = Hist('TH2F', sample=sa, var=['true_phi', 'gen_phi'], binning=(24, -3.142, 3.142, 24, -3.142, 3.142), sel=X.get('$' + sel), wt='wt_pu*wt_def')
                    hists[year]['true_reco_phi'][sa][sel] = Hist('TH2F', sample=sa, var=['true_phi', 'reco_phi'], binning=(24, -3.142, 3.142, 24, -3.142, 3.142), sel=X.get('$' + sel), wt='wt_pu*wt_def')
                    # hists[year]['sphi_response'][sa][sel] = Hist('TH2F', sample=sa, var=['gen_sphi', 'reco_sphi'], binning=(6, -1.571, 1.571, 6, -1.571, 1.571), sel=X.get('$' + sel), wt='wt_pu*wt_def')
    # X.get('$photon_pix', printlevel=1)
    for year in years:
        MultiDraw(hists[year], sample_files[year], 'WGDataAnalysis', mt_cores=4, mt_thresh=1E5)

        for path, hname, obj in hists[year].ListObjects():
            name = obj.sample
            if name is not 'data_obs':
                tgt_lumi = sample_cfgs[year][samples[year]['data_obs']]['lumi']
                events = sample_cfgs[year][samples[year][name]]['events']
                xsec = sample_cfgs[year][samples[year][name]]['xsec']
                scale = tgt_lumi * xsec / events
                obj.Scale(scale)

        for sel in ['photon_pix_m', 'photon_pix_e', 'eft_m', 'eft_e']:
            for sa in ['WG']:
                NormTH2InColumns(hists[year]['sphi_response'][sa][sel])
                NormTH2InColumns(hists[year]['sphi_xy_response'][sa][sel])
                NormTH2InColumns(hists[year]['sphi_puppi_response'][sa][sel])


    fout = ROOT.TFile('output_notebooks.root', 'RECREATE')
    NodeToTDir(fout, hists)
    fout.Close()
else:
    fin = ROOT.TFile('output_notebooks.root')
    TDirToNode(fin, node=hists)

x_titles = {
    'l0_pt': ['Gen. lepton p_{T}', 'GeV'],
    'p0_pt': ['Gen. photon p_{T}', 'GeV'],
    'lhe_p0_pt': ['LHE photon p_{T}', 'GeV'],
    'lhe_p0_pt_wide': ['LHE photon p_{T}', 'GeV'],
}
plotcfg = DEFAULT_CFG
plotcfg.update({
    'type': 'multihist',
    'legend_pos': [0.55, 0.66, 0.90, 0.91],
    'ratio_pad_frac': 0.33,
    'top_title_right': '41.5 fb^{-1} (13 TeV, 2017)',
    'logy': True,
    'logy_min': 1E-4
    })


# Draw efficiencies for the stitched sample
for year, var in itertools.product(years, ['l0_pt', 'p0_pt']):
    hist_dict = {}
    for opath, objname, obj in hists[year][var]['WG'].ListObjects(depth=0):
        hist_dict[objname] = obj

    # Draw N/N-1 efficiencies
    MakeMultiHistPlot('efficiencies_m_%s_%s' % (year, var),
                      outdir=outdir,
                      hists=hist_dict,
                      cfg=UpdateDict(plotcfg, {
                          'x_title': x_titles[var],
                          'ratio_y_title': '#varepsilon(N/N-1)',
                          'top_title_right': '%.1f fb^{-1} (13 TeV, %s)' % (lumi_fb[year], year),
                          'ratio_y_range': [0.55, 1.05],
                          }),
                      layout=[
                          {'name': 'baseline_m', 'legend': 'Baseline'},
                          {'name': 'muon_id', 'legend': ' + Muon ID'},
                          {'name': 'muon_iso', 'legend': ' + Muon Iso'},
                          {'name': 'muon_trg', 'legend': ' + Muon Trg(24/27) OR (50)'},
                          {'name': 'muon_trg_2', 'legend': ' + Muon Trg(24/27)'},
                          {'name': 'photon_id_m', 'legend': ' + Photon ID'},
                          {'name': 'photon_pix_m', 'legend': ' + Photon pix veto'},
                      ],
                      ratios=[
                          {'num': 'muon_id', 'den': 'baseline_m', 'type': 'binomial'},
                          {'num': 'muon_iso', 'den': 'muon_id', 'type': 'binomial'},
                          {'num': 'muon_trg', 'den': 'muon_iso', 'type': 'binomial'},
                          {'num': 'muon_trg_2', 'den': 'muon_trg', 'type': 'binomial'},
                          {'num': 'photon_id_m', 'den': 'muon_trg_2', 'type': 'binomial'},
                          {'num': 'photon_pix_m', 'den': 'photon_id_m', 'type': 'binomial'},
                      ]
    )
    MakeMultiHistPlot('efficiencies_e_%s_%s' % (year, var),
                      outdir=outdir,
                      hists=hist_dict,
                      cfg=UpdateDict(plotcfg, {
                          'x_title': x_titles[var],
                          'ratio_y_title': '#varepsilon(N/N-1)',
                          'top_title_right': '%.1f fb^{-1} (13 TeV, %s)' % (lumi_fb[year], year),
                          'ratio_y_range': [0.45, 1.05],
                          }),
                      layout=[
                          {'name': 'baseline_e', 'legend': 'Baseline'},
                          {'name': 'elec_id', 'legend': ' + Electron ID'},
                          # {'name': 'elec_match', 'legend': ' + Electron matched'},
                          # {'name': 'elec_iso', 'legend': ' + Electron Iso'},
                          {'name': 'elec_trg', 'legend': ' + Electron Trg'},
                          # {'name': 'elec_trg_2', 'legend': ' + Electron Trg'},
                          {'name': 'photon_id_e', 'legend': ' + Photon ID'},
                          # {'name': 'photon_match_e', 'legend': ' + Photon matched'},
                          {'name': 'photon_pix_e', 'legend': ' + Photon pix veto'},
                      ],
                      ratios=[
                          {'num': 'elec_id', 'den': 'baseline_e', 'type': 'binomial'},
                          # {'num': 'elec_match', 'den': 'elec_id', 'type': 'binomial'},
                          # {'num': 'elec_iso', 'den': 'elec_id', 'type': 'binomial'},
                          {'num': 'elec_trg', 'den': 'elec_id', 'type': 'binomial'},
                          # {'num': 'elec_trg_2', 'den': 'elec_trg', 'type': 'binomial'},
                          {'num': 'photon_id_e', 'den': 'elec_trg', 'type': 'binomial'},
                          # {'num': 'photon_match_e', 'den': 'photon_id_e', 'type': 'binomial'},
                          {'num': 'photon_pix_e', 'den': 'photon_id_e', 'type': 'binomial'},
                      ]
    )

    # Draw N/baseline efficiencies
    MakeMultiHistPlot('efficiencies_m_cumulative_%s_%s' % (year, var),
                      outdir=outdir,
                      hists=hist_dict,
                      cfg=UpdateDict(plotcfg, {
                          'x_title': x_titles[var],
                          'ratio_y_title': '#varepsilon(N/baseline)',
                          'top_title_right': '%.1f fb^{-1} (13 TeV, %s)' % (lumi_fb[year], year),
                          'ratio_y_range': [0.35, 1.05],
                      }),
                      layout=[
                          {'name': 'baseline_m', 'legend': 'Baseline'},
                          {'name': 'muon_id', 'legend': ' + Muon ID'},
                          {'name': 'muon_iso', 'legend': ' + Muon Iso'},
                          {'name': 'muon_trg', 'legend': ' + Muon Trg(24/27) OR (50)'},
                          {'name': 'muon_trg_2', 'legend': ' + Muon Trg(24/27)'},
                          {'name': 'photon_id_m', 'legend': ' + Photon ID'},
                          {'name': 'photon_pix_m', 'legend': ' + Photon pix veto'},
                      ],
                      ratios=[
                          {'num': 'muon_id', 'den': 'baseline_m', 'type': 'binomial'},
                          {'num': 'muon_iso', 'den': 'baseline_m', 'type': 'binomial'},
                          {'num': 'muon_trg', 'den': 'baseline_m', 'type': 'binomial'},
                          {'num': 'muon_trg_2', 'den': 'baseline_m', 'type': 'binomial'},
                          {'num': 'photon_id_m', 'den': 'baseline_m', 'type': 'binomial'},
                          {'num': 'photon_pix_m', 'den': 'baseline_m', 'type': 'binomial'},
                      ]
    )
    MakeMultiHistPlot('efficiencies_e_cumulative_%s_%s' % (year, var),
                      outdir=outdir,
                      hists=hist_dict,
                      cfg=UpdateDict(plotcfg, {
                          'x_title': x_titles[var],
                          'ratio_y_title': '#varepsilon(N/baseline)',
                          'top_title_right': '%.1f fb^{-1} (13 TeV, %s)' % (lumi_fb[year], year),
                          'ratio_y_range': [0.35, 1.05],
                      }),
                      layout=[
                          {'name': 'baseline_e', 'legend': 'Baseline'},
                          {'name': 'elec_id', 'legend': ' + Electron ID'},
                          # {'name': 'elec_iso', 'legend': ' + Electron Iso'},
                          {'name': 'elec_trg', 'legend': ' + Electron Trg'},
                          # {'name': 'elec_trg_2', 'legend': ' + Electron Trg'},
                          {'name': 'photon_id_e', 'legend': ' + Photon ID'},
                          {'name': 'photon_pix_e', 'legend': ' + Photon pix veto'},
                      ],
                      ratios=[
                          {'num': 'elec_id', 'den': 'baseline_e', 'type': 'binomial'},
                          # {'num': 'elec_iso', 'den': 'baseline_e', 'type': 'binomial'},
                          {'num': 'elec_trg', 'den': 'baseline_e', 'type': 'binomial'},
                          # {'num': 'elec_trg_2', 'den': 'baseline_e', 'type': 'binomial'},
                          {'num': 'photon_id_e', 'den': 'baseline_e', 'type': 'binomial'},
                          {'num': 'photon_pix_e', 'den': 'baseline_e', 'type': 'binomial'},
                      ]
    )


# Draw LO vs stitched vs NLO per year
# The LO cross sections should have been adjusted to give the same normalisation as
# the NLO under the baseline selection above
for year, var in itertools.product(years, ['lhe_p0_pt', 'lhe_p0_pt_wide', 'p0_pt', 'l0_pt']):
    hist_dict = {}
    for hname in ['WG', 'WG-inc', 'WG-NLO']:
        hist_dict[hname] = hists[year][var][hname]['baseline_m']

    cfg = deepcopy(plotcfg)
    cfg.update({
        'x_title': x_titles[var],
        'ratio_y_range': [0.65, 1.35],
        'ratio_y_title': 'Ratio to LO'
        })

    MakeMultiHistPlot('lhe_%s_%s' % (year, var),
                      outdir=outdir,
                      hists=hist_dict,
                      cfg=UpdateDict(plotcfg, {
                          'x_title': x_titles[var],
                          'ratio_y_range': [0.65, 1.35],
                          'ratio_y_title': 'Ratio to LO',
                          'legend_show_yields': True
                      }),
                      layout=[
                          {'name': 'WG', 'legend': 'LO stitched sample'},
                          {'name': 'WG-inc', 'legend': 'LO inclusive sample'},
                          {'name': 'WG-NLO', 'legend': 'NLO inclusive sample'}
                      ],
                      ratios=[
                          {'num': 'WG-inc', 'den': 'WG'},
                          {'num': 'WG-NLO', 'den': 'WG'}
                      ])

    MakeMultiHistPlot('lhe_normed_%s_%s' % (year, var),
                      outdir=outdir,
                      hists=hist_dict,
                      cfg=UpdateDict(plotcfg, {
                          'x_title': x_titles[var],
                          'ratio_y_range': [0.65, 1.35],
                          'ratio_y_title': 'Ratio to LO',
                          'legend_show_yields': True,
                          'norm_to': 1.0,
                          'top_title_right': '%.1f fb^{-1} (13 TeV, %s)' % (lumi_fb[year], year)
                      }),
                      layout=[
                          {'name': 'WG', 'legend': 'LO stitched sample'},
                          {'name': 'WG-NLO', 'legend': 'NLO inclusive sample'}
                      ],
                      ratios=[
                          {'num': 'WG-NLO', 'den': 'WG'}
                      ])
