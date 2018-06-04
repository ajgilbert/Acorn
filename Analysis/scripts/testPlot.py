import ROOT
import time
import json
from pprint import pprint
from collections import defaultdict

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)


class Node:
    def __init__(self):
        self.d = {}

    def __getitem__(self, key):
        if key in self.d:
            return self.d[key]
        else:
            self.d[key] = Node()
            return self.d[key]

    def __setitem__(self, key, val):
        self.d[key] = val

    def ListNodes(self, depth=9999999, res=None, dirname=''):
        if res is None:
            res = []
        pre = '/' if dirname is not '' else ''
        res.append((dirname, self))
        if depth > 0:
            for key, val in self.d.iteritems():
                if isinstance(val, Node):
                    val.ListNodes(depth - 1, res, dirname + pre + key)
        return res

    def ListObjects(self, depth=9999999, res=None, dirname=''):
        if res is None:
            res = []
        nodelist = self.ListNodes(depth=depth, dirname=dirname)
        for path, node in nodelist:
            for key, val in node.d.iteritems():
                if not isinstance(val, Node):
                    res.append((path, key, val))
        return res


def WriteToTFile(obj, filedir, path, name):
    filedir.cd()
    startdir = ROOT.gDirectory.GetPath()
    as_vec = path.split('/')
    if len(as_vec) >= 1:
        for i in xrange(0, len(as_vec)):
            if not ROOT.gDirectory.GetDirectory(as_vec[i]):
                ROOT.gDirectory.mkdir(as_vec[i])
            ROOT.gDirectory.cd(as_vec[i])
    if not ROOT.gDirectory.FindKey(name):
        obj.SetName(name)
        ROOT.gDirectory.WriteTObject(obj, name)
    ROOT.gDirectory.cd(startdir)

def Define(classname, args, sample, var, wt='1.0', sel='1'):
    res = getattr(ROOT, classname)("", "", *args)
    res.classname = classname
    res.sample = sample
    res.var = var
    res.wt = wt
    res.sel = sel
    return res


def GenTreeCode(tree, objlist):
    fname = 'RunTree%i' % GenTreeCode.counter
    # static variable to keep ensure we can give
    # each function a unique name
    GenTreeCode.counter += 1
    res = []
    res.append('void %s(TTree *tree, TObjArray *hists) {' % fname)
    res.append('  TTreeReader reader(tree);')
    branches = tree.GetListOfBranches()
    binfos = []
    for i in range(branches.GetEntriesFast()):
        branch = branches.UncheckedAt(i)
        bname = branch.GetName()
        btype = branch.GetListOfLeaves().UncheckedAt(0).GetTypeName()
        binfos.append((bname, btype))
        res.append('  TTreeReaderValue<%s> %s_(reader, "%s");' % (btype, bname, bname))
    res.append('  while (reader.Next()) {')
    for bname, btype in binfos:
        res.append('    %s %s = *%s_;' % (btype, bname, bname))
    for i, obj in enumerate(objlist):
        res.append('    if (%s) static_cast<%s*>(hists->UncheckedAt(%i))->Fill(%s, %s);' % (obj.sel, obj.classname, i, ', '.join(obj.var), obj.wt))
    res.append('  }')
    res.append('}')
    print '\n'.join(res)
    start = time.time()
    ROOT.gInterpreter.Declare('\n'.join(res))
    end = time.time()
    print '>> Getting JIT-y with it for function %s in %.2g seconds' % (fname, (end - start))
    return fname


GenTreeCode.counter = 0

tname = 'DiMuonAnalysis'

# samples = {
#     'DY': 'output/020618/wgamma_2016_v1/Main/DYJetsToLL_M-50-madgraphMLM.root',
#     'data_obs': 'output/020618/wgamma_2016_v1/Main/SingleMuon.root',
#     'TT': 'output/020618/wgamma_2016_v1/Main/TT-powheg.root',
#     'WG': 'output/020618/wgamma_2016_v1/Main/WGToLNuG-madgraphMLM.root'
# }

# remap = {
#     'DY': 'DYJetsToLL_M-50-madgraphMLM',
#     'data_obs': 'SingleMuon',
#     'TT': 'TT-powheg',
#     'WG': 'WGToLNuG-madgraphMLM'
# }

samples = {
    'data_obs': 'output/020618/wgamma_2018_v1/Main/SingleMuon.root',
    'DY': 'output/020618/wgamma_2018_v1/Main/DYJetsToLL_M-50-madgraphMLM.root',
    'TT': 'output/020618/wgamma_2018_v1/Main/TTTo2L2Nu-powheg.root'
}

remap = {
    'DY': 'DYJetsToLL_M-50-madgraphMLM',
    'data_obs': 'SingleMuon',
    'TT': 'TTTo2L2Nu-powheg',
    'WG': 'WGToLNuG-madgraphMLM'
}

hists = Node()

for ptcut in range(30, 100, 5):
    for sample in samples:
        hists['m_ll']['pt_%s' % ptcut][sample] = Define('TH1F', (60, 60., 120.), sample, ['m_ll'], sel='pt_1 > %s' % ptcut)


drawtree = defaultdict(list)
for path, name, obj in hists.ListObjects():
    drawtree[obj.sample].append(obj)

pprint(drawtree)
with open('input/cfg_wgamma_2018_v1.json') as jsonfile:
    cfg = json.load(jsonfile)
    sample_cfg = cfg['samples']

for sample, drawobjs in drawtree.iteritems():
    f = ROOT.TFile(samples[sample])
    sample_cfg[remap[sample]]['events'] = f.Get('counters').GetBinContent(2)
    t = f.Get(tname)
    fname = GenTreeCode(t, drawobjs)
    histarr = ROOT.TObjArray(len(drawobjs))
    for i, h in enumerate(drawobjs):
        histarr.AddAt(h, i)
    start = time.time()
    getattr(ROOT, fname)(t, histarr)
    end = time.time()
    runtime = end - start
    print '>> Evaluated sample %s with %i events in %.2g seconds (%.2g evts/sec)' % (sample, t.GetEntriesFast(), runtime, t.GetEntriesFast() / runtime)
    f.Close()


for path, name, obj in hists.ListObjects():
    if name is not 'data_obs':
        tgt_lumi = sample_cfg[remap['data_obs']]['lumi']
        events = sample_cfg[remap[name]]['events']
        xsec = sample_cfg[remap[name]]['xsec']
        scale = tgt_lumi * xsec / events
        obj.Scale(scale)

fout = ROOT.TFile('output_2018.root', 'RECREATE')

for path, name, obj in hists.ListObjects():
    WriteToTFile(obj, fout, path, name)

fout.Close()
