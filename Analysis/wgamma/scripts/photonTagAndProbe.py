import ROOT
import argparse
from array import array
from Acorn.Analysis.analysis import *

ROOT.RooWorkspace.imp = getattr(ROOT.RooWorkspace, 'import')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(0)

parser = argparse.ArgumentParser()
parser.add_argument('-y', '--year', default='2016')
args = parser.parse_args()

bin_cfgs = [
    {
        'name': 'WorstIso_pt_bins_inc_eta',
        'var': 'm_ll',
        'binning': (50, 75, 125),
        'tag': 't_trg && t_id && t_rand',
        'probe': 'p_id',
        'binvar_x': 'p_pt',
        'bins_x': [30., 40., 50., 60., 100., 200., 1000.],
        'binvar_y': 'abs(p_eta)',
        'bins_y': [0, 2.5]
    },
    {
        'name': 'WorstIso_pt_eta_bins',
        'var': 'm_ll',
        'binning': (50, 75, 125),
        'tag': 't_trg && t_id && t_rand',
        'probe': 'p_id',
        'binvar_x': 'p_pt',
        'bins_x': [30., 40., 50., 60., 100., 200., 1000.],
        'binvar_y': 'abs(p_eta)',
        'bins_y': [0, 1.0, 1.56, 2.1, 2.5]
    }
]

drawlist = []
# andable = set()

for cfg in bin_cfgs:
    cfg['hist'] = ROOT.TH2D(cfg['name'], cfg['name'],
                            len(cfg['bins_x'])-1, array('d', cfg['bins_x']),
                            len(cfg['bins_y'])-1, array('d', cfg['bins_y']))
    hist = cfg['hist']
    hist.GetXaxis().SetTitle(cfg['binvar_x'])
    hist.GetYaxis().SetTitle(cfg['binvar_y'])

    cfg['bins'] = []

    for i in xrange(1, hist.GetNbinsX()+1):
        for j in xrange(1, hist.GetNbinsY()+1):
            cfg['bins'].append('%s>=%g && %s<%g && %s>=%g && %s<%g' % (
                cfg['binvar_x'], hist.GetXaxis().GetBinLowEdge(i),
                cfg['binvar_x'], hist.GetXaxis().GetBinUpEdge(i),
                cfg['binvar_y'], hist.GetYaxis().GetBinLowEdge(j),
                cfg['binvar_y'], hist.GetYaxis().GetBinUpEdge(j),
                ))
            # andable.add('%s>=%g' % (cfg['binvar_x'], hist.GetXaxis().GetBinLowEdge(i)))
            # andable.add('%s<%g' % (cfg['binvar_x'], hist.GetXaxis().GetBinUpEdge(i)))
            # andable.add('%s>=%g' % (cfg['binvar_y'], hist.GetYaxis().GetBinLowEdge(j)))
            # andable.add('%s<%g' % (cfg['binvar_y'], hist.GetYaxis().GetBinUpEdge(j)))

    for b in cfg['bins']:
        drawlist.append((cfg['var'], '((%s) && !(%s) && (%s)) * wt' % (b, cfg['probe'], cfg['tag'])))
        drawlist.append((cfg['var'], '((%s) && (%s) && (%s)) * wt' % (b, cfg['probe'], cfg['tag'])))
        # andable.add(cfg['probe'])
        # andable.add(cfg['tag'])


remaps = {
    "2016": {
        'DYJetsToLL': 'DYJetsToLL_M-50-amcatnloFXFX',
        # 'data_obs_m': 'SingleMuon',
        'Data': 'SingleElectron',
    },
    "2017": {
        'DYJetsToLL': 'DYJetsToLL_M-50-amcatnloFXFX',
        # 'data_obs_m': 'SingleMuon',
        'Data': 'SingleElectron',
    },
    "2018": {
        'DYJetsToLL': 'DYJetsToLL_M-50-amcatnloFXFX',
        # 'data_obs_m': 'SingleMuon',
        'Data': 'EGamma',
    }
}

remap = remaps[args.year]
prefix = 'root://eoscms.cern.ch//store/user/agilbert/ANv5-200511-full/wgamma_%s_v5/PhotonTP_' % args.year

samples = {}
for sa in remap:
    samples[sa] = (prefix + remap[sa] + '.root')


hists = Node()


# trees = {
#     'DYJetsToLL': analysis.TTreeEvaluator('ZeeTP', 'output/HTT2016Studies_July19/DYJetsToLL.root'),
#     'Data': analysis.TTreeEvaluator('ZeeTP', 'output/HTT2016Studies_July19/SingleElectron.root')
# }


# sys.exit(0)

for sample in remap:
    outfile = ROOT.TFile('PhotonTP_%s_%s.root' % (args.year, sample), 'RECREATE')

    hists = Node()
    for cfg in bin_cfgs:

        for b in cfg['bins']:
            hists[cfg['name']]['%s:fail' % b] = Hist('TH1F', sample=sample, var=[cfg['var']], binning=cfg['binning'], sel='((%s) && !(%s) && (%s))' % (b, cfg['probe'], cfg['tag']), wt='wt_def')
            hists[cfg['name']]['%s:pass' % b] = Hist('TH1F', sample=sample, var=[cfg['var']], binning=cfg['binning'], sel='((%s) && (%s) && (%s))' % (b, cfg['probe'], cfg['tag']), wt='wt_def')
            # drawlist.append((cfg['var'], '((%s) && !(%s) && (%s)) * wt' % (b, cfg['probe'], cfg['tag'])))
            # drawlist.append((cfg['var'], '((%s) && (%s) && (%s)) * wt' % (b, cfg['probe'], cfg['tag'])))

    MultiDraw(hists, samples, 'WGTagAndProbe', mt_cores=4)

    # hists = trees[sample].Draw(drawlist, compiled=True)

    i = 0
    NodeToTDir(outfile, hists)
    for cfg in bin_cfgs:
        wsp = ROOT.RooWorkspace('wsp_'+cfg['name'], '')
        var = wsp.factory('m_ll[100,75,125]')

        # outfile.cd()
        # outfile.mkdir(cfg['name'])
        # ROOT.gDirectory.cd(cfg['name'])

        for b in cfg['bins']:
            # hists[2*i].SetName(b+':fail')
            # hists[2*i+1].SetName(b+':pass')
            # hists[2*i].Write()
            # hists[2*i+1].Write()
            # hists[cfg['name']]['%s:fail' % b].Write()
            # hists[cfg['name']]['%s:pass' % b].Write()
            dat = wsp.imp(ROOT.RooDataHist(b, '', ROOT.RooArgList(var),
                          ROOT.RooFit.Index(wsp.factory('cat[fail,pass]')),
                          ROOT.RooFit.Import('fail', hists[cfg['name']]['%s:fail' % b]),
                          ROOT.RooFit.Import('pass', hists[cfg['name']]['%s:pass' % b]))
                          )
            i += 1
        outfile.cd()
        wsp.Write()
        cfg['hist'].Write()
        wsp.Delete()

    outfile.Close()
