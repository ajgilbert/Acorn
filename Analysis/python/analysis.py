import ROOT
import time
# import json
from array import array
from pprint import pprint
from collections import defaultdict, OrderedDict

ROOT.TH1.AddDirectory(False)

class Node:
    def __init__(self):
        self.d = OrderedDict()

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

    def ForEach(self, action, depth=9999999):
        objs = self.ListObjects(depth)
        for path, name, obj in objs:
            action(obj)


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


def Hist(classname, binning, sample, var, wt='1.0', sel='1'):
    # If a list is supplied interpret as variable binning
    if isinstance(binning, list):
        binargs = (len(binning) - 1, array('d', binning))
        res = getattr(ROOT, classname)("", "", *binargs)
    else:
        res = getattr(ROOT, classname)("", "", *binning)
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
    # print '\n'.join(res)
    start = time.time()
    ROOT.gInterpreter.Declare('\n'.join(res))
    end = time.time()
    print '>> JIT compiled function %s in %.2g seconds' % (fname, (end - start))
    return fname


GenTreeCode.counter = 0


def MultiDraw(node, sample_to_file_dict, tree_name):
    drawtree = defaultdict(list)
    codegen_time = 0
    draw_time = 0
    for path, name, obj in node.ListObjects():
        drawtree[obj.sample].append(obj)
    # pprint(drawtree)
    for sample, drawobjs in drawtree.iteritems():
        f = ROOT.TFile(sample_to_file_dict[sample])
        t = f.Get(tree_name)
        start = time.time()
        fname = GenTreeCode(t, drawobjs)
        end = time.time()
        codegen_time += (end - start)
        histarr = ROOT.TObjArray(len(drawobjs))
        for i, h in enumerate(drawobjs):
            histarr.AddAt(h, i)
        start = time.time()
        getattr(ROOT, fname)(t, histarr)
        end = time.time()
        runtime = end - start
        draw_time += runtime
        print '>> Evaluated sample %s with %i events in %.2g seconds (%.2g evts/sec)' % (sample, t.GetEntriesFast(), runtime, t.GetEntriesFast() / runtime)
        f.Close()
    print '>> Total code generation time: %.2g seconds' % codegen_time
    print '>> Total evaluation time: %.2g seconds' % draw_time


def NormaliseTo(hist, val):
    hist.Scale(1. / hist.Integral())
