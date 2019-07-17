import subprocess
import json
import glob
import argparse
from Acorn.Analysis.analysis import *

parser = argparse.ArgumentParser()
parser.add_argument('--config', default='fid_pt_binned')
parser.add_argument('--steps', default='')
parser.add_argument('--years', default='2016,2017,2018')
args = parser.parse_args()

def call(args):
    print subprocess.list2cmdline(args)
    subprocess.call(args)


configs = {
    "fid_pt_binned": {
        'pt_bins': '[30,40,50,60,80,100,200,500]',
        # 'pt_bins': '[30,40,50,60,80,100,120,160,200,250,300,500]',
        'phi_var': '0.5',
        'phi_var_label': '1',
        'phi_var_obs': '0.5',
        'phi_bins': '(1,0.,1.)',
        'phi_bins_obs': '(1,0.,1.)',
        'task_name': 'fid_region'
    },
    "phi_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi)',
        'phi_var_label': '|#phi|',
        'phi_var_obs': 'abs(reco_phi)',
        'phi_bins': '(5,0.,math.pi)',
        'phi_bins_obs': '(5,0.,math.pi)',
        'task_name': 'eft_region'
    },
    "phi_f_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)',
        'task_name': 'eft_region'
    },
    "puppi_phi_f_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_puppi_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)',
        'task_name': 'eft_region'
    },
    "pt_binned": {
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
    for sample, sample_type in [
        ('/home/files/190411-gen/wgamma_2016_v2/wg_gen_WGToMuNuG-EFT-madgraphMLM-stitched.root', 'LO'),
        ('/home/files/190411-gen/wgamma_2016_v3/wg_gen_WGToMuNuG_01J_5f_EFT-stitched.root', 'NLO')
    ]:

        plot_dir = '/eos/user/a/agilbert/www/wgamma/eft_scaling/%s/%s' % (label, sample_type)
        call(['mkdir', '-p', plot_dir])

        standard_cuts = [
            '--g_pt', '150', '--l_pt', '80.0',
            '--l_eta', '2.4', '--g_eta', '2.5', '--dr', '3.0',
            '--sample', sample,
            '--plot-dir', plot_dir
        ]

        for chg in ['+1', '-1']:
            call(['python', 'wgamma/scripts/makeEFTScaling.py',
                  '--draw-x', 'g_pt', pt_bins, 'p_{T}^{#gamma}',
                  '--draw-y', phi_var, phi_bins, phi_var_label,
                  '--n_pt', '80.0',
                  '--charge', chg, '--o', '%s_%s' % (label, sample_type), '--label', 'main'] + standard_cuts)
            call(['python', 'wgamma/scripts/makeEFTScaling.py',
                  '--draw-x', 'g_pt', pt_bins, 'p_{T}^{#gamma}',
                  '--draw-y', phi_var, phi_bins, phi_var_label,
                  '--n_pt', '40.0', '--n_pt_max', '80.0',
                  '--charge', chg, '--o', '%s_%s' % (label, sample_type), '--label', 'met1'] + standard_cuts)

    call(['gallery.py', '/eos/user/a/agilbert/www/wgamma/eft_scaling'])


if 'makeHists' in steps:
    testplot_args = {
        'pt_bins': pt_bins,
        'phi_var': phi_var,
        'phi_var_obs': phi_var_obs,
        'phi_bins': phi_bins,
        'phi_bins_obs': phi_bins_obs
    }
    print json.dumps(testplot_args)
    for yr in years:
        indir = 'root://eoscms.cern.ch//store/cmst3/user/agilbert/190630-full/wgamma_%s_v4/WGamma_' % yr
        # call(['python', 'wgamma/scripts/makeHists.py', '--task', config['task_name'],
        #       '--indir', indir,
        #       '--year', yr, '--extra-cfg', json.dumps(testplot_args), '--label', label])
        do_systs = [
          ('MetJesLo_', '_CMS_scale_met_jesDown'),
          ('MetJesHi_', '_CMS_scale_met_jesUp'),
          # ('PScaleLo_', '_CMS_scale_pDown'),
          # ('PScaleHi_', '_CMS_scale_pUp'),
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
                  ])

if 'T2W' in steps:
    infiles = []
    for region in ['main', 'met1']:
        for sgn in ['p', 'n']:
            infiles.append('%s_LO_%s_%s.root:%s:%s' % (using_label, region, sgn, region, sgn))

    call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s/*.txt' % label)) +
         ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
          '--PO', 'ptBins=%i' % n_pt_bins,
          '--PO', 'phiBins=%i' % n_phi_bins,
          '--PO', 'type=eft', '--PO',
          'files=%s' % ','.join(infiles), '--channel-masks', '-o', 'combined_%s.root' % label])

    # if label == 'fid_pt_binned':
    #     call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s/*.txt' % label)) +
    #          ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
    #           '--PO', 'ptBins=%i' % n_pt_bins,
    #           '--PO', 'phiBins=%i' % n_phi_bins,
    #           '--PO', 'type=pt_diff', '-o', 'combined_pt_diff_%s.root' % label])
    # else:
    #     call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s/*.txt' % label)) +
    #          ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
    #           '--PO', 'ptBins=%i' % n_pt_bins,
    #           '--PO', 'phiBins=%i' % n_phi_bins,
    #           '--PO', 'type=pt_phi_diff', '-o', 'combined_pt_phi_diff_%s.root' % label])

if 'limitsVsPtMax' in steps:
    for bsm_label, bsm_setting in [('withBSM', 1), ('noBSM', 0)]:
        for i in range(n_pt_bins):
            setpars = ['c3w=0', 'lumiscale=1', 'withBSM=%i' % bsm_setting]
            setpars = ','.join(setpars + ['mask_%i=1' % X for X in range(i + 1, n_pt_bins)])
            call(['combine', '-M', 'AsymptoticLimits', '-t', '-1', '-n', '.%s.%s' % (label, bsm_label),
                  '-m', '%i' % i, '--setParameters', setpars, 'combined_%s.root' % label,
                  '--setParameterRanges', 'c3w=0,0.1'])
        call(['combineTool.py', '-M', 'CollectLimits', '-o', 'limits_%s_%s.json' % (label, bsm_label)] +
              glob.glob('higgsCombine.%s.%s.*.root' % (label, bsm_label)))


allPOIs = []
for i in range(n_pt_bins):
    if label == 'fid_pt_binned':
        allPOIs.append('r_x_%i' % i)
    else:
        for j in range(n_phi_bins):
            allPOIs.append('r_p_%i_%i' % (i, j))
            allPOIs.append('r_n_%i_%i' % (i, j))

if 'xsec2D' in steps:
    initPOIs = ','.join([('%s=1' % X) for X in allPOIs])
    genStr = 'P;n;;' + ';'.join(['%s,%s' % (X, X) for X in allPOIs])

    if label == 'fid_pt_binned':
        wsp = 'pt_diff'
        rangePOIs = ':'.join([('%s=0.5,1.5' % X) for X in allPOIs])
    else:
        wsp = 'pt_phi_diff'
        rangePOIs = ':'.join([('%s=0.5,1.5' % X) for X in allPOIs])

    call(['combineTool.py', '-M', 'MultiDimFit', 'combined_%s_%s.root' % (wsp, label), '-t', '-1',
          '--algo', 'grid', '--setParameters', initPOIs, '--setParameterRanges', rangePOIs, '--floatOtherPOIs', '1',
          '--generate', genStr, '-n', '.%s' % label, '--points', '30', '--alignEdges', '1', '--parallel', '4'])

if 'xsec2DPlot' in steps:

    call(['rm', '%s.json' % label])

    for POI in allPOIs:
        call(['python', 'wgamma/scripts/plot1DScan.py',
              '--main', 'higgsCombine.%s.%s.MultiDimFit.mH120.root' % (label, POI), '--POI', POI,
              '--model', 'xsec2D', '--output', 'scan_%s_%s' % (label, POI),
              '--json', '%s.json' % label, '--chop', '20'])


"""
python wgamma/scripts/theoryLimitPlot.py limits_phi_binned_noBSM.json limits_pt_binned_noBSM.json:exp0:'MarkerSize=0,LineWidth=2,LineColor=4,Title="No #phi binning"' --show exp --limit-on "C_{3W} (TeV^{-2})" --x-title "Maximum p_{T}^{#gamma} (GeV)"
"""
