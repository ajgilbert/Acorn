import subprocess
import json
import glob
from Acorn.Analysis.analysis import *


def call(args):
    print subprocess.list2cmdline(args)
    subprocess.call(args)


configs = {
    "phi_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi)',
        'phi_var_label': '|#phi|',
        'phi_var_obs': 'abs(reco_phi)',
        'phi_bins': '(5,0.,math.pi)',
        'phi_bins_obs': '(5,0.,math.pi)'
    },
    "phi_f_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)'
    },
    "puppi_phi_f_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_puppi_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)',
        'phi_bins_obs': '(3,0.,math.pi/2.)'
    },
    "pt_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi)',
        'phi_var_label': '|#phi|',
        'phi_var_obs': 'p0_pt',
        'phi_bins': '(5,0.,math.pi)',
        'phi_bins_obs': '(5,0.,math.pi)',
        'using_label': 'phi_binned'
    },
}

label = 'phi_f_binned'
using_label = label if 'using_label' not in configs[label] else configs[label]['using_label']
pt_bins = configs[label]['pt_bins']
phi_var = configs[label]['phi_var']
phi_var_label = configs[label]['phi_var_label']
phi_var_obs = configs[label]['phi_var_obs']
phi_bins = configs[label]['phi_bins']

n_pt_bins = len(BinEdgesFromStr(pt_bins)) - 1
n_phi_bins = len(BinEdgesFromStr(phi_bins)) - 1

steps = [
    # 'makeEFTScaling',
    'makeHists',
    # 'setupDatacards',
    # 'T2W',
    # 'limitsVsPtMax'
    # 'xsec2D',
    # 'xsec2DPlot'
]

years = [
    '2016',
    # '2017',
    # '2018'
]

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
        'phi_bins': phi_bins
    }
    print json.dumps(testplot_args)
    for yr in years:
        call(['python', 'wgamma/scripts/makeHists.py', '--task', 'eft_region',
              '--indir', '/home/files/190412-full/wgamma_%s_v3/WGamma_' % yr,
              '--year', yr, '--extra-cfg', json.dumps(testplot_args), '--label', label])


if 'setupDatacards' in steps:
    for yr in years:
        for chn in ['e', 'm']:
            call(['python', 'wgamma/scripts/setupDatacards.py', '--var', phi_var_obs,
                  '--output', 'output/cards/%s' % label, '--label', using_label,
                  '--year', yr, '--channel', chn, '--type', 'eft',
                  '--pt-bins', '%i' % n_pt_bins,
                  '--phi-bins', '%i' % n_phi_bins
                  ])

if 'T2W' in steps:
    infiles = []
    for region in ['main', 'met1']:
        for sgn in ['p', 'n']:
            infiles.append('%s_LO_%s_%s.root:%s:%s' % (using_label, region, sgn, region, sgn))

    # call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s/*.txt' % label)) +
    #      ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
    #       '--PO', 'ptBins=%i' % n_pt_bins,
    #       '--PO', 'phiBins=%i' % n_phi_bins,
    #       '--PO', 'type=eft', '--PO',
    #       'files=%s' % ','.join(infiles), '--channel-masks', '-o', 'combined_%s.root' % label])

    # call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s/*.txt' % label)) +
    #      ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
    #       '--PO', 'ptBins=%i' % n_pt_bins,
    #       '--PO', 'phiBins=%i' % n_phi_bins,
    #       '--PO', 'type=pt_diff', '-o', 'combined_pt_diff_%s.root' % label])

    call(['combineTool.py', '-M', 'T2W', '-i'] + list(glob.glob('output/cards/%s/*.txt' % label)) +
         ['--cc', '-P', 'Acorn.Analysis.WGPhysicsModel:wgModel',
          '--PO', 'ptBins=%i' % n_pt_bins,
          '--PO', 'phiBins=%i' % n_phi_bins,
          '--PO', 'type=pt_phi_diff', '-o', 'combined_pt_phi_diff_%s.root' % label])

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


if 'xsec2D' in steps:
    allPOIs = []
    for i in range(n_pt_bins):
        for j in range(n_phi_bins):
            allPOIs.append('r_p_%i_%i' % (i, j))
            allPOIs.append('r_n_%i_%i' % (i, j))
    initPOIs = ','.join([('%s=1' % X) for X in allPOIs])

    genStr = 'P;n;;' + ';'.join(['%s,%s' % (X, X) for X in allPOIs])

    call(['combineTool.py', '-M', 'MultiDimFit', 'combined_pt_phi_diff_%s.root' % label, '-t', '-1',
          '--algo', 'grid', '--setParameters', initPOIs, '--floatOtherPOIs', '1',
          '--generate', genStr, '-n', '.%s' % label, '--points', '30', '--alignEdges', '1', '--parallel', '4'])

if 'xsec2DPlot' in steps:
    allPOIs = []
    for i in range(n_pt_bins):
        for j in range(n_phi_bins):
            allPOIs.append('r_p_%i_%i' % (i, j))
            allPOIs.append('r_n_%i_%i' % (i, j))

    call(['rm', '%s.json' % label])

    for POI in allPOIs:
        call(['python', 'wgamma/scripts/plot1DScan.py',
              '--main', 'higgsCombine.%s.%s.MultiDimFit.mH120.root' % (label, POI), '--POI', POI,
              '--model', 'xsec2D', '--output', 'scan_%s_%s' % (label, POI),
              '--json', '%s.json' % label, '--chop', '20'])


"""
python wgamma/scripts/theoryLimitPlot.py limits_phi_binned_noBSM.json limits_pt_binned_noBSM.json:exp0:'MarkerSize=0,LineWidth=2,LineColor=4,Title="No #phi binning"' --show exp --limit-on "C_{3W} (TeV^{-2})" --x-title "Maximum p_{T}^{#gamma} (GeV)"
"""
