import ROOT
import json
# import sys
from pprint import pprint
from collections import defaultdict
import argparse
from Acorn.Analysis.analysis import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)


parser = argparse.ArgumentParser()

parser.add_argument('--task', default='eft_region', choices=['eft_region', 'photon_fakes'])
parser.add_argument('--indir', default='output/130818-reduced/wgamma_2016_v2/WGamma/')

args = parser.parse_args()


tname = 'WGDataAnalysis'
prefix = args.indir

remap = {
    'DY': 'DYJetsToLL_M-50-madgraphMLM',
    'data_obs': 'SingleMuon',
    'TT': 'TT-powheg',
    'WG': 'WGToLNuG-madgraphMLM-stitched',
    'W': 'WJetsToLNu-madgraphMLM',
    'VVTo2L2Nu': 'VVTo2L2Nu-amcatnloFXFX',
    'WWTo1L1Nu2Q': 'WWTo1L1Nu2Q-amcatnloFXFX',
    'WZTo1L1Nu2Q': 'WZTo1L1Nu2Q-amcatnloFXFX',
    'WZTo1L3Nu': 'WZTo1L3Nu-amcatnloFXFX',
    'WZTo2L2Q': 'WZTo2L2Q-amcatnloFXFX',
    'WZTo3LNu': 'WZTo3LNu-amcatnloFXFX',
    'ZZTo2L2Q': 'ZZTo2L2Q-amcatnloFXFX',
    'ZZTo4L': 'ZZTo4L-amcatnloFXFX'
}

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')


hists = Node()
do_cats = []


X = SelectionManager()

if args.task == 'eft_region':
    pt_bins_min = [150, 210, 300, 420]
    pt_bins_max = [210, 300, 420, 600]
    phi_bins_min = [0.00, 0.63, 1.26, 1.89, 2.52]
    phi_bins_max = [0.63, 1.26, 1.89, 2.52, 3.15]

    X.Set('baseline', sel='m0_trg && n_m==1 && m0_iso<0.15 && n_p==1 && p0_medium && !p0_haspix && m0_pt>80 && met>80 && p0_pt>150 && m0p0_dr>3.0', wt='wt_def*wt_pu*wt_m0*wt_trg_m0*wt_p0')
    do_cats.extend(['baseline'])

    X.Set('p_gen_ooa', sel='gen_m0_q==+1 && !(gen_m0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_m0p0_dr>3.0)')
    X.Set('n_gen_ooa', sel='gen_m0_q==-1 && !(gen_m0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_m0p0_dr>3.0)')

    X.Set('p_gen_acc', sel='gen_m0_q==+1 && gen_m0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_m0p0_dr>3.0')
    X.Set('n_gen_acc', sel='gen_m0_q==-1 && gen_m0_pt>80 && gen_met>80 && gen_p0_pt>150 && gen_m0p0_dr>3.0')

    for i in range(len(pt_bins_min)):
        # Reconstructed categories
        X.Derive('p_%i' % i, 'baseline', sel='m0_q==+1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i]))
        X.Derive('n_%i' % i, 'baseline', sel='m0_q==-1 && p0_pt>=%f && p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i]))
        do_cats.extend(['p_%i' % i, 'n_%i' % i])

        # Gen level selections
        X.Derive('p_gen_%i' % (i), 'p_gen_acc', sel='gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i]))
        X.Derive('n_gen_%i' % (i), 'n_gen_acc', sel='gen_p0_pt>=%f && gen_p0_pt<%f' % (pt_bins_min[i], pt_bins_max[i]))

        for j in range(len(phi_bins_min)):
            X.Derive('p_gen_%i_%i' % (i, j), 'p_gen_acc', sel='gen_p0_pt>=%f && gen_p0_pt<%f && abs(gen_reco_phi) >= %f && abs(gen_reco_phi) < %f' % (pt_bins_min[i], pt_bins_max[i], phi_bins_min[j], phi_bins_max[j]))
            X.Derive('n_gen_%i_%i' % (i, j), 'n_gen_acc', sel='gen_p0_pt>=%f && gen_p0_pt<%f && abs(gen_reco_phi) >= %f && abs(gen_reco_phi) < %f' % (pt_bins_min[i], pt_bins_max[i], phi_bins_min[j], phi_bins_max[j]))

    drawvars = [
        ('p0_pt', pt_bins_min + [pt_bins_max[-1]]),
        ('abs(reco_phi)', (5, 0, 3.15)),
        ('abs(gen_phi)', (5, 0, 3.15)),
    ]

    for sel in do_cats:
        for var, binning in drawvars:
            for sample in samples:
                hists[sel][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel=X.sel('$' + sel), wt=X.wt('$' + sel))
            for P in ['W', 'DY', 'TT']:
                hists[sel][var]['%s_R' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.sel('$' + sel + ' && p0_isprompt'), wt=X.wt('$' + sel))
                hists[sel][var]['%s_F' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.sel('$' + sel + ' && !p0_isprompt'), wt=X.wt('$' + sel))
            hists[sel][var]['WG_p_ooa'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $p_gen_ooa'), wt=X.wt('$' + sel))
            hists[sel][var]['WG_n_ooa'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $n_gen_ooa'), wt=X.wt('$' + sel))
            hists[sel][var]['WG_p_acc'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $p_gen_acc'), wt=X.wt('$' + sel))
            hists[sel][var]['WG_n_acc'] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $n_gen_acc'), wt=X.wt('$' + sel))
            for i in range(len(pt_bins_min)):
                hists[sel][var]['WG_p_%i' % (i)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $p_gen_%i' % (i)), wt=X.wt('$' + sel))
                hists[sel][var]['WG_n_%i' % (i)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $n_gen_%i' % (i)), wt=X.wt('$' + sel))
                for j in range(len(phi_bins_min)):
                    hists[sel][var]['WG_p_%i_%i' % (i, j)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $pho_p_gen_%i_%i' % (i, j)), wt=X.wt('$' + sel))
                    hists[sel][var]['WG_n_%i_%i' % (i, j)] = Hist('TH1D', sample='WG', var=[var], binning=binning, sel=X.sel('$' + sel + ' && $pho_n_gen_%i_%i' % (i, j)), wt=X.wt('$' + sel))


if args.task == 'photon_fakes':
    X.Set('baseline', sel='m0_trg && n_m==1 && m0_iso<0.15 && m0met_mt>60 && n_p==1 && m0p0_dr>0.7 && p0_medium_noch && !p0_haspix', wt='wt_def*wt_pu*wt_m0*wt_trg_m0*wt_p0')
    X.Derive('barrel', 'baseline', sel='abs(p0_eta) < 1.4442')
    X.Derive('endcap', 'baseline', sel='abs(p0_eta) > 1.4442')
    do_cats.extend(['baseline', 'barrel', 'endcap'])
    sdb = {
        'barrel': '0.01022',
        'endcap': '0.03001'
    }
    for S in ['barrel', 'endcap']:
        X.Derive('%s_iso_l' % S, S, sel='p0_chiso > 2.5 && p0_chiso < 7')
        X.Derive('%s_sig_l' % S, S, sel='p0_sigma > %s' % sdb[S])
        X.Derive('%s_iso_l_sig_t' % S, S, sel='p0_chiso > 2.5 && p0_chiso < 7 && p0_sigma < %s' % sdb[S])
        X.Derive('%s_iso_l_sig_l' % S, S, sel='p0_chiso > 2.5 && p0_chiso < 7 && p0_sigma > %s' % sdb[S])
        X.Derive('%s_iso_t_sig_l' % S, S, sel='p0_chiso < 0.441 && p0_sigma > %s' % sdb[S])
        do_cats.extend(['%s_iso_l' % S, '%s_sig_l' % S, '%s_iso_l_sig_t' % S, '%s_iso_l_sig_l' % S, '%s_iso_t_sig_l' % S])

    drawvars = [
        ('m0met_mt', (30, 0., 200.)),
        ('m0_pt', (40, 0., 150.)),
        ('m0_eta', (20, -3.0, 3.0)),
        ('m0_iso', (40, 0, 2.0)),
        ('m1_pt', (40, 0., 150.)),
        ('m1_eta', (20, -3.0, 3.0)),
        ('m0m1_M', (40, 60, 120)),
        ('m0m1_dr', (20, 0., 5.)),
        ('met', (20, 0., 200.)),
        ('p0_pt', (20, 0, 200.)),
        ('p0_eta', (20, -3.0, 3.0)),
        ('m0p0_dr', (20, 0., 5.)),
        ('m0p0_M', (40, 60, 120)),
        ('p0_chiso', (40, 0, 20.0)),
        ('p0_neiso', (40, 0, 20.0)),
        ('p0_phiso', (40, 0, 20.0)),
        ('p0_hovere', (20, 0., 0.5)),
        ('p0_sigma', (40, 0., 0.050)),
        ('p0_haspix', (2, -0.5, 1.5))
    ]

    for sel in do_cats:
        for var, binning in drawvars:
            for sample in samples:
                hists[sel][var][sample] = Hist('TH1D', sample=sample, var=[var], binning=binning, sel=X.sel('$'+sel), wt=X.wt('$'+sel))
                for P in ['W', 'DY', 'TT']:
                    hists[sel][var]['%s_R' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.sel('$' + sel + ' && p0_isprompt'), wt=X.wt('$' + sel))
                    hists[sel][var]['%s_F' % P] = Hist('TH1D', sample=P, var=[var], binning=binning, sel=X.sel('$' + sel + ' && !p0_isprompt'), wt=X.wt('$' + sel))


MultiDraw(hists, samples, tname)

with open('input/cfg_wgamma_2016_v2.json') as jsonfile:
    cfg = json.load(jsonfile)
    sample_cfg = cfg['samples']

for sample in samples:
    f = ROOT.TFile(samples[sample])
    sample_cfg[remap[sample]]['events'] = f.Get('counters').GetBinContent(2)
    f.Close()

for path, hname, obj in hists.ListObjects():
    name = obj.sample
    if name is not 'data_obs':
        tgt_lumi = sample_cfg[remap['data_obs']]['lumi']
        events = sample_cfg[remap[name]]['events']
        xsec = sample_cfg[remap[name]]['xsec']
        scale = tgt_lumi * xsec / events
        obj.Scale(scale)

fout = ROOT.TFile('output_2016_%s.root' % args.task, 'RECREATE')

for path, name, obj in hists.ListObjects():
    WriteToTFile(obj, fout, path, name)

fout.Close()
