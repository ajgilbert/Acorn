import subprocess
import json
import glob
import argparse
import os
from Acorn.Analysis.analysis import *
from Acorn.Analysis.jobs import Jobs

job_mgr = Jobs()
parser = argparse.ArgumentParser()
job_mgr.attach_job_args(parser)
parser.add_argument('--config', default='fid_pt_binned')
parser.add_argument('--steps', default='')
parser.add_argument('--years', default='2016,2017,2018')
args = parser.parse_args()
job_mgr.set_args(args)


def call(args, noSub=False):
    print subprocess.list2cmdline(args)
    if noSub:
        subprocess.call(args)
    else:
        job_mgr.job_queue.append(subprocess.list2cmdline(args))


configs = {
    "fid_pt_binned": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[30,50,70,100,150,200,300,500,800,1200]',
        # 'pt_bins': '[30,40,50,60,80,100,120,160,200,250,300,500]',
        'phi_var': '0.5',
        'phi_var_label': '1',
        'phi_var_obs': '0.5',
        'phi_bins': '(1,0.,1.)',
        'phi_bins_obs': '(1,0.,1.)',
        'task_name': 'fid_region'
    },
    "fid_p0_eta_binned": {
        'x_var': 'gen_p0_eta',
        'x_var_obs': 'p0_eta',
        'pt_bins': '[-2.5,-1.9,-1.5,-1.1,-0.75,-0.45,-0.15,0.15,0.45,0.75,1.1,1.5,1.9,2.5]',
        'phi_var': '0.5',
        'phi_var_label': '1',
        'phi_var_obs': '0.5',
        'phi_bins': '(1,0.,1.)',
        'phi_bins_obs': '(1,0.,1.)',
        'task_name': 'fid_region'
    },
    "fid_mt_cluster_binned": {
        'x_var': 'gen_mt_cluster',
        'x_var_obs': 'mt_cluster',
        'pt_bins': '[0,100,150,200,250,300,400,500,650,800,1000,1200,1700,2500]',
        'phi_var': '0.5',
        'phi_var_label': '1',
        'phi_var_obs': '0.5',
        'phi_bins': '(1,0.,1.)',
        'phi_bins_obs': '(1,0.,1.)',
        'task_name': 'fid_region'
    },
    "fid_l0p0_dr_binned": {
        'x_var': 'gen_l0p0_dr',
        'x_var_obs': 'l0p0_dr',
        'pt_bins': '[0.7,1.0,1.3,1.6,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.0,4.5,5.0]',
        'phi_var': '0.5',
        'phi_var_label': '1',
        'phi_var_obs': '0.5',
        'phi_bins': '(1,0.,1.)',
        'phi_bins_obs': '(1,0.,1.)',
        'task_name': 'fid_region'
    },
    "inclusive_xs": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[30,10000]',
        # 'pt_bins': '[30,40,50,60,80,100,120,160,200,250,300,500]',
        'phi_var': '0.5',
        'phi_var_label': '1',
        'phi_var_obs': '0.5',
        'phi_bins': '(1,0.,1.)',
        'phi_bins_obs': '(1,0.,1.)',
        'task_name': 'fid_region'
    },
    "phi_binned": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi)',
        'phi_var_label': '|#phi|',
        'phi_var_obs': 'abs(reco_phi)',
        'phi_bins': '(5,0.,math.pi)',
        'phi_bins_obs': '(5,0.,math.pi)',
        'task_name': 'eft_region'
    },
    "phi_f_binned": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)',
        'task_name': 'eft_region'
    },
    "puppi_phi_f_binned": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[150,200,300,500,800,1200]',
        'phi_var': 'abs(gen_true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_puppi_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)',
        'task_name': 'eft_region',
        'jet_veto': False,
        'split_charge': True
    },
    "puppi_phi_f_binned_nochg": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[150,200,300,500,800,1200]',
        'phi_var': 'abs(gen_true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_puppi_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)',
        'task_name': 'eft_region',
        'jet_veto': False,
        'split_charge': False
    },
    "puppi_phi_f_binned_jetveto": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[150,200,300,500,800,1200]',
        'phi_var': 'abs(gen_true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_puppi_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)',
        'task_name': 'eft_region',
        'jet_veto': True
    },
    "puppi_phi_f_binned_nobin": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[150,200,300,500,800,1200]',
        'phi_var': 'abs(true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'p0_pt',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)',
        'task_name': 'eft_region',
        'using_label': 'puppi_phi_f_binned'
    },
    "pt_binned": {
        'x_var': 'gen_p0_pt',
        'x_var_obs': 'p0_pt',
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi)',
        'phi_var_label': '|#phi|',
        'phi_var_obs': 'p0_pt',
        'phi_bins': '(5,0.,math.pi)',
        'phi_bins_obs': '(5,0.,math.pi)',
        'using_label': 'phi_binned',
        'task_name': 'eft_region'
    },
}

label = args.config
config = configs[label]

using_label = label if 'using_label' not in config else config['using_label']
x_var = config['x_var']
x_var_obs = config['x_var_obs']
pt_bins = config['pt_bins']
phi_var = config['phi_var']
phi_var_label = config['phi_var_label']
phi_var_obs = config['phi_var_obs']
phi_bins = config['phi_bins']
phi_bins_obs = config['phi_bins_obs']

n_pt_bins = len(BinEdgesFromStr(pt_bins)) - 1
n_phi_bins = len(BinEdgesFromStr(phi_bins)) - 1

# steps = [
#     # 'makeEFTScaling',
#     'makeHists',
#     # 'setupDatacards',
#     # 'T2W',
#     # 'limitsVsPtMax'
#     # 'xsec2D',
#     # 'xsec2DPlot'
# ]
steps = args.steps.split(',')

years = args.years.split(',')
# years = [
#     # '2016',
#     # '2017',
#     '2018'
# ]

if 'makeEFTScaling' in steps:
    outdir = 'eft_scaling_new_2018_st'
    for sample, sample_type in [
        # ('/home/files/190411-gen/wgamma_2016_v2/wg_gen_WGToMuNuG-EFT-madgraphMLM-stitched.root', 'LO'),
        # ('/home/files/190411-gen/wgamma_2016_v3/wg_gen_WGToMuNuG_01J_5f_EFT-stitched.root', 'NLO')
        # ('/eos/cms/store/user/agilbert/ANv5-200430-gen/wgamma_2016_v5/wg_gen_WGToMuNuG-EFT-madgraphMLM-stitched.root', 'LO'),
        # ('/eos/cms/store/user/agilbert/ANv5-200430-gen/wgamma_2016_v5/wg_gen_WGToMuNuG_01J_5f_EFT-stitched.root', 'NLO')
        # ('/eos/cms/store/user/agilbert/ANv5-200430-gen/wgamma_2018_v5/wg_gen_WGToLNuG-madgraphMLM-stitched.root', 'LO'),
        # ('/eos/cms/store/user/agilbert/ANv5-200622-gen/wgamma_2018_v5/wg_gen_WGToLNuG-madgraphMLM-stitched.root', 'LO'),
        # ('/eos/cms/store/user/agilbert/ANv5-200430-gen/wgamma_2018_v5/wg_gen_WGToLNuG-amcatnloFXFX-stitched.root', 'NLO')
        ('/eos/cms/store/user/agilbert/ANv5-200624-gen/wgamma_2016_v5/wg_gen_WGToMuNuG_01J_5f_EFT-stitched.root', 'NLO')
    ]:

        plot_dir = '/eos/user/a/agilbert/www/wgamma/%s/%s/%s' % (outdir, label, sample_type)
        call(['mkdir', '-p', plot_dir])

        standard_cuts = [
            '--g_pt', '150', '--l_pt', '80.0',
            '--l_eta', '2.4', '--g_eta', '2.5', '--dr', '3.0',
            '--sample', sample,
            '--plot-dir', plot_dir
        ]

        standard_cuts += ['--pre-wt', 'wt', '--pre-var', 'gen']

        charges = ['0']
        if 'split_charge' in config and config['split_charge']:
            charges = ['+1', '-1']

        for chg in charges:
            call(['python', 'wgamma/scripts/makeEFTScaling.py',
                  '--draw-x', 'gen_p0_pt', pt_bins, 'p_{T}^{#gamma}',
                  '--draw-y', phi_var, phi_bins, phi_var_label,
                  '--n_pt', '80.0',
                  '--charge', chg, '--o', '%s_%s' % (label, sample_type), '--label', 'main'] + standard_cuts, noSub=True)
            call(['python', 'wgamma/scripts/makeEFTScaling.py',
                  '--draw-x', 'gen_p0_pt', pt_bins, 'p_{T}^{#gamma}',
                  '--draw-y', phi_var, phi_bins, phi_var_label,
                  '--n_pt', '40.0', '--n_pt_max', '80.0',
                  '--charge', chg, '--o', '%s_%s' % (label, sample_type), '--label', 'met1'] + standard_cuts, noSub=True)

    call(['gallery.py', '/eos/user/a/agilbert/www/wgamma/%s' % outdir], noSub=True)

if 'makeHists' in steps:
    testplot_args = {
        'x_var': x_var,
        'x_var_obs': x_var_obs,
        'pt_bins': pt_bins,
        'phi_var': phi_var,
        'phi_var_obs': phi_var_obs,
        'phi_bins': phi_bins,
        'phi_bins_obs': phi_bins_obs,
        'jet_veto': config['jet_veto'],
        'split_charge': config['split_charge']
    }
    print json.dumps(testplot_args)
    for yr in years:
        indir = 'root://eoscms.cern.ch//store/user/agilbert/ANv5-200623-full/wgamma_%s_v5/WGamma_' % yr
        call(['python', 'wgamma/scripts/makeHists.py', '--task', config['task_name'],
              '--indir', indir,
              '--year', yr, '--extra-cfg', json.dumps(testplot_args), '--label', label])
        do_systs = [
          ('MetJesLo_', '_CMS_scale_met_jesDown'),
          ('MetJesHi_', '_CMS_scale_met_jesUp'),
          ('MetUncLo_', '_CMS_scale_met_unclusteredDown'),
          ('MetUncHi_', '_CMS_scale_met_unclusteredUp'),
          ('PScaleLo_', '_CMS_scale_pDown'),
          ('PScaleHi_', '_CMS_scale_pUp'),
        ]
        for syst_file, syst_name in do_systs:
            call(['python', 'wgamma/scripts/makeHists.py', '--task', config['task_name'],
                  '--indir-data', indir,
                  '--indir', indir + syst_file, '--syst', syst_name,
                  '--year', yr, '--extra-cfg', json.dumps(testplot_args), '--label', label + syst_name])


if 'setupDatacards' in steps:
    for yr in years:
        for chn in ['e', 'm']:
            call(['python', 'wgamma/scripts/setupDatacards.py', '--var', phi_var_obs,
                  '--output', 'output/cards/%s' % label, '--label', using_label,
                  '--year', yr, '--channel', chn, '--type', config['task_name'],
                  '--pt-bins', '%i' % n_pt_bins,
                  '--phi-bins', '%i' % n_phi_bins
                  ], noSub=True)

if 'T2W' in steps:
    infiles = []
    for region in ['main', 'met1']:
        charges = ['x']
        if config['split_charge']:
            charges = ['p', 'n']
        for sgn in [charges]:
            # infiles.append('%s_NLO_%s_%s.root:%s:%s' % (using_label, region, sgn, region, sgn))
            if label == 'puppi_phi_f_binned_jetveto':
                infiles.append('CMS_2020_PAS_SMP_20_005_eft_%s_%s_photon_pt_phi_jveto.json:%s:%s' % (region, sgn, region, sgn))
            else:
                infiles.append('CMS_2020_PAS_SMP_20_005_eft_%s_%s_photon_pt_phi.json:%s:%s' % (region, sgn, region, sgn))

    if label == 'fid_pt_binned' or label == 'inclusive_xs':
        call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s' % label)) +
             ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
              '--PO', 'ptBins=%i' % n_pt_bins,
              '--PO', 'phiBins=%i' % n_phi_bins,
              '--PO', 'type=pt_diff', '-o', '%s/combined_pt_diff_%s.root' % (os.getcwd(), label)], noSub=True)
    else:
        call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s' % label)) +
             ['--cc', 'combined_EFT.txt', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
              '--PO', 'ptBins=%i' % n_pt_bins,
              '--PO', 'phiBins=%i' % n_phi_bins,
              '--PO', 'type=eft', '--PO',
              'files=%s' % ','.join(infiles), '--channel-masks', '-o', '%s/combined_%s.root' % (os.getcwd(), label)], noSub=True)

        call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s' % label)) +
             ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
              '--PO', 'ptBins=%i' % n_pt_bins,
              '--PO', 'phiBins=%i' % n_phi_bins,
              '--PO', 'type=pt_phi_diff', '-o', '%s/combined_pt_phi_diff_%s.root' % (os.getcwd(), label)], noSub=False)

if 'limitsVsPtMax' in steps:
    for bsm_label, bsm_setting, fit_range in [
        ('withBSM', 1, 1.5),
        ('noBSM', 0, 6.0)
    ]:
        for i in range(n_pt_bins):
            setpars = ['c3w=0', 'lumiscale=1', 'withBSM=%i' % bsm_setting]
            setpars = ','.join(setpars + ['mask_%i=1' % X for X in range(i + 1, n_pt_bins)])
            par_range = fit_range * (1. - float(i) / float(n_pt_bins))
            # call(['combine', '-M', 'AsymptoticLimits', '-t', '-1', '-n', '.%s.%s' % (label, bsm_label),
            #       '-m', '%i' % i, '--setParameters', setpars, 'combined_%s.root' % label,
            #       '--setParameterRanges', 'c3w=0,0.1'])
            call(['combine', '-M', 'MultiDimFit', '-t', '-1', '-n', '.%s.%s' % (label, bsm_label),
                  '--algo', 'grid', '-m', '%i' % i, '--setParameters', setpars, 'combined_%s.root' % label,
                  '--setParameterRanges', 'c3w=%f,%f' % (-1. * par_range, par_range), '--points', '40', '--alignEdges', '1',
                  '--cminDefaultMinimizerStrategy', '0', '--X-rtd', 'OPTIMIZE_BOUNDS=0'])
            # call(['combine', '-M', 'MultiDimFit', '-t', '1', '-s', '1', '-n', '.%s.%s.toy' % (label, bsm_label),
            #       '--algo', 'grid', '-m', '%i' % i, '--setParameters', setpars, 'combined_%s.root' % label,
            #       '--setParameterRanges', 'c3w=%f,%f' % (-1. * par_range, par_range), '--points', '40', '--alignEdges', '1',
            #       '--cminDefaultMinimizerStrategy', '0', '--X-rtd', 'MINIMIZER_analytic', '--X-rtd', 'OPTIMIZE_BOUNDS=0'])


allPOIs = []
for i in range(n_pt_bins):
    if label == 'fid_pt_binned' or label == 'inclusive_xs':
        allPOIs.append('r_x_%i' % i)
    else:
        for j in range(n_phi_bins):
            allPOIs.append('r_p_%i_%i' % (i, j))
            allPOIs.append('r_n_%i_%i' % (i, j))

if 'xsec2D' in steps:
    initPOIs = ','.join([('%s=1' % X) for X in allPOIs])
    genStr = 'P;n;;' + ';'.join(['%s,%s' % (X, X) for X in allPOIs])
    if label == 'fid_pt_binned' or label == 'inclusive_xs':
        wsp = 'pt_diff'
        rangePOIs = ':'.join([('%s=0.5,1.5' % X) for X in allPOIs])
        grpStr = 'freezeNuisanceGroups;n;;!,,nominal;th,fr.th;th,expt,,fr.expt;all,fr.all;all,autoMCStats,,fr.mcstats'
    else:
        wsp = 'pt_phi_diff'
        rangePOIs = ':'.join([('%s=0,10' % X) for X in allPOIs])
        # Only do a more basic breakdown for the 2D
        grpStr = 'freezeNuisanceGroups;n;;sig_inc,,nominal;all,autoMCStats,,fr.mcstats'

    call(['combineTool.py', '-M', 'MultiDimFit', 'combined_%s_%s.root' % (wsp, label), '-t', '-1',
          '--algo', 'grid', '--setParameters', initPOIs, '--setParameterRanges', rangePOIs, '--floatOtherPOIs', '1',
          '--generate', genStr, grpStr,
          '-n', '.%s' % label, '--points', '30', '--alignEdges', '1', '--parallel', '4',
          '--cminDefaultMinimizerStrategy', '0', '--X-rtd', 'MINIMIZER_analytic', '--X-rtd', 'OPTIMIZE_BOUNDS=0'])
    # call(['combineTool.py', '-M', 'MultiDimFit', 'combined_%s_%s.root' % (wsp, label), '-t', '-1',
    #       '--algo', 'none', '--setParameters', initPOIs, '--setParameterRanges', rangePOIs, '--floatOtherPOIs', '1',
    #       '-n', '.%s.corr' % label, '--parallel', '4', '--saveFitResult', '--robustHesse', '0', '-v', '3',
    #       '--cminDefaultMinimizerStrategy', '0', '--X-rtd', 'MINIMIZER_analytic', '--X-rtd', 'OPTIMIZE_BOUNDS=0'])

if 'xsec2DPlot' in steps:

    call(['rm', '%s.json' % label])
    call(['rm', '%s_breakdown.json' % label])

    for POI in allPOIs:
        call(['python', 'wgamma/scripts/plot1DScan.py',
            '--main', 'higgsCombine.%s.%s.nominal.MultiDimFit.mH120.root' % (label, POI), '--POI', POI,
            '--model', 'xsec2D', '--output', 'scan_%s_%s' % (label, POI),
            '--json', '%s.json' % label, '--chop', '20', '--translate', 'wgamma/scripts/pois.json',
            '--others',
            "higgsCombine.%s.%s.fr.mcstats.MultiDimFit.mH120.root:Freeze all:2" % (label, POI),
            '--breakdown', 'Syst,Stat'])

        if label == 'fid_pt_binned' or label == 'inclusive_xs':
            call(['python', 'wgamma/scripts/plot1DScan.py',
                  '--main', 'higgsCombine.%s.%s.nominal.MultiDimFit.mH120.root' % (label, POI), '--POI', POI,
                  '--model', 'xsec2D', '--output', 'scan_%s_%s_breakdown' % (label, POI),
                  '--json', '%s_breakdown.json' % label, '--chop', '20', '--translate', 'wgamma/scripts/pois.json',
                  '--others',
                  "higgsCombine.%s.%s.fr.th.MultiDimFit.mH120.root:Freeze Th.:2"  % (label, POI),
                  "higgsCombine.%s.%s.fr.expt.MultiDimFit.mH120.root:Freeze Th. + Expt.:4"  % (label, POI),
                  "higgsCombine.%s.%s.fr.all.MultiDimFit.mH120.root:Freeze all:8"  % (label, POI),
                  "higgsCombine.%s.%s.fr.mcstats.MultiDimFit.mH120.root:Freeze all + MCStats:12" % (label, POI),
                  '--breakdown', 'Th,Expt,Lumi,MCStats,Stat'])

if 'C3WPlot' in steps:
    for bsm_label in ['withBSM', 'noBSM']:
        for i in range(n_pt_bins):
            call(['python', 'wgamma/scripts/plot1DScan.py', '--limit-cls',
                '--main', 'higgsCombine.%s.%s.MultiDimFit.mH%i.root' % (label, bsm_label, i), '--POI', 'c3w',
                '--model', 'c3w_bin_%i' % i, '--output', 'scan_%s_%s_%i' % (label, bsm_label, i),
                '--json', '%s_%s.json' % (label, bsm_label), '--chop', '30', '--remove-near-min', '0.8'])
            # call(['python', 'wgamma/scripts/plot1DScan.py',
            #     '--main', 'higgsCombine.%s.%s.toy.MultiDimFit.mH%i.1.root' % (label, bsm_label, i), '--POI', 'c3w',
            #     '--model', 'c3w_bin_%i' % i, '--output', 'scan_toy_%s_%s_%i' % (label, bsm_label, i),
            #     '--json', '%s_%s_toy.json' % (label, bsm_label), '--chop', '30', '--remove-near-min', '0.8'])
"""
python wgamma/scripts/theoryLimitPlot.py limits_phi_binned_noBSM.json limits_pt_binned_noBSM.json:exp0:'MarkerSize=0,LineWidth=2,LineColor=4,Title="No #phi binning"' --show exp --limit-on "C_{3W} (TeV^{-2})" --x-title "Maximum p_{T}^{#gamma} (GeV)"
"""

job_mgr.flush_queue()
