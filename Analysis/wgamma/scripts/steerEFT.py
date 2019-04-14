import subprocess
import json
from Acorn.Analysis.analysis import *


def call(args):
    print subprocess.list2cmdline(args)
    subprocess.call(args)


label = 'phi_binned'
pt_bins = '[150,210,300,420,600,850,1200]'
phi_var = 'abs(true_phi)'
phi_var_obs = 'abs(reco_phi)'
phi_bins = '(5,0.,math.pi)'

# label = 'phi_f_binned'
# pt_bins = '[150,210,300,420,600,850,1200]'
# phi_var = 'abs(true_phi_f)'
# phi_bins = '(3,0.,math.pi/2.)'

steps = [
    # 'genPlot2',
    'testPlotWG',
    # 'setupDatacards'
]

years = [
    '2016',
    # '2017',
    # '2018'
]


if 'genPlot2' in steps:
    standard_cuts = [
        '--g_pt', '150', '--l_pt', '80.0',
        '--l_eta', '2.4', '--g_eta', '2.5', '--dr', '3.0'
    ]
    for chg in ['+1', '-1']:
        call(['python', 'wgamma/scripts/genPlot2.py',
              '--draw-x', 'g_pt', pt_bins, 'pt',
              '--draw-y', phi_var, phi_bins, 'phi',
              '--n_pt', '80.0',
              '--charge', chg, '--o', 'main_%s' % label] + standard_cuts)
        call(['python', 'wgamma/scripts/genPlot2.py',
              '--draw-x', 'g_pt', pt_bins, 'pt',
              '--draw-y', phi_var, phi_bins, 'phi',
              '--n_pt', '40.0', '--n_pt_max', '80.0',
              '--charge', chg, '--o', 'met1_%s' % label] + standard_cuts)


if 'testPlotWG' in steps:
    testplot_args = {
        'pt_bins': pt_bins,
        'phi_var': phi_var,
        'phi_var_obs': phi_var_obs,
        'phi_bins': phi_bins
    }
    print json.dumps(testplot_args)
    for yr in years:
        call(['python', 'wgamma/scripts/testPlotWG.py', '--task', 'eft_region',
              '--indir', '/home/files/190406-full/wgamma_%s_v3/WGamma_' % yr,
              '--year', yr, '--extra-cfg', json.dumps(testplot_args), '--label', label])


if 'setupDatacards' in steps:
    for yr in years:
        for chn in ['e', 'm']:
            call(['python', 'wgamma/scripts/setupDatacards.py', '--var', phi_var_obs,
                  '--output', 'output/cards/%s' % label, '--label', label,
                  '--year', yr, '--channel', chn, '--type', 'eft',
                  '--pt-bins', '%i' % (len(BinEdgesFromStr(pt_bins)) - 1),
                  '--phi-bins', '%i' % (len(BinEdgesFromStr(phi_bins)) - 1)
                  ])
