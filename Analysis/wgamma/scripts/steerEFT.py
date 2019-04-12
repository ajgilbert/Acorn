import subprocess
import math


def BinningArgStr(bins):
    if isinstance(bins, list):
        return '[%s]' % ','.join([str(x) for x in bins])
    if isinstance(bins, tuple):
        return '(%i,%f,%f)' % bins
        # width = (float(bins[2]) - float(bins[1])) / float(bins[0])
        # binedges = []
        # for i in range(bins[0]):
        #     binedges.append(float(bins[1]) + float(i) * width)


pt_bins = [150, 210, 300, 420, 600, 850, 1200]
phi_var = 'abs(true_phi_f)'
phi_bins = (3, 0., math.pi / 2.)

print BinningArgStr(pt_bins)
print BinningArgStr(phi_bins)

standard_cuts = [
    '--g_pt', '150', '--l_pt', '80.0',
    '--l_eta', '2.4', '--g_eta', '2.5', '--dr', '3.0'
]
for chg in ['+1', '-1']:
    subprocess.call(['python', 'wgamma/scripts/genPlot2.py',
                    '--draw-x', 'g_pt', BinningArgStr(pt_bins), 'pt',
                    '--draw-y', phi_var, BinningArgStr(phi_bins), 'phi',
                    '--n_pt', '80.0',
                    '--charge', chg, '--o', 'main'] + standard_cuts)
    subprocess.call(['python', 'wgamma/scripts/genPlot2.py',
                    '--draw-x', 'g_pt', BinningArgStr(pt_bins), 'pt',
                    '--draw-y', phi_var, BinningArgStr(phi_bins), 'phi',
                    '--n_pt', '40.0', '--n_pt_max', '80.0',
                    '--charge', chg, '--o', 'met1'] + standard_cuts)
