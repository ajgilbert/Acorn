import subprocess
import json
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
        'phi_bins': '(5,0.,math.pi)'
    },
    "phi_f_binned": {
        'pt_bins': '[150,210,300,420,600,850,1200]',
        'phi_var': 'abs(true_phi_f)',
        'phi_var_label': '|#phi_{f}|',
        'phi_var_obs': 'abs(reco_phi_f)',
        'phi_bins': '(3,0.,math.pi/2.)'
    }
}

label = 'phi_binned'
pt_bins = configs[label]['pt_bins']
phi_var = configs[label]['phi_var']
phi_var_label = configs[label]['phi_var_label']
phi_var_obs = configs[label]['phi_var_obs']
phi_bins = configs[label]['phi_bins']

steps = [
    'makeEFTScaling',
    'makeHists',
    'setupDatacards'
]

years = [
    '2016',
    '2017',
    '2018'
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
              '--indir', '/home/files/190412-full/wgamma_%s_v3/WGamma_' % yr,
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

"""

combineTool.py -M T2W -i output/cards/phi_binned/p_*.txt --cc -P Acorn.Analysis.WGPhysicsModel:wgModel \
--PO ptBins=6 --PO phiBins=5 --PO type=eft --PO files=phi_binned_LO_main_p.root:main,phi_binned_LO_met1_p.root:met1 --channel-masks -o combined_phi_binned.root

combineTool.py -M T2W -i output/cards/phi_f_binned/p_*.txt --cc -P Acorn.Analysis.WGPhysicsModel:wgModel \

--PO ptBins=6 --PO phiBins=3 --PO type=eft --PO files=main_phi_f_binned_w_p.root:main,met1_phi_f_binned_w_p.root:met1 --channel-masks -o combined_phi_f_binned.root




"""
