import ROOT
import time
import re
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


class SelectionManager:
    def __init__(self):
        self.storage = OrderedDict()
        self.varregex = re.compile('\$[_a-zA-Z][_a-zA-Z0-9]*')
        self.idx = 0

    def Set(self, label, sel='1', wt='1.'):
        self.storage[label] = (sel, wt)

    def Derive(self, label, base, sel='1', wt=None):
        if wt is None:
            new_wt = self.storage[base][1]
        else:
            new_wt = '$%s * (%s)' % (base, wt)
        self.Set(label, sel='$%s && (%s)' % (base, sel), wt=new_wt)

    def _getreplacement(self, matchgroup, idx):
        label = matchgroup.group(0)[1:]
        return '(' + self.storage[label][idx] + ')'

    def _doreplacement(self, pattern, expr, idx=0):
        newexpr = pattern.sub(lambda x: self._getreplacement(x, idx), expr)
        # print '%s --> %s' % (expr, newexpr)
        if newexpr == expr:
            return newexpr
        else:
            return self._doreplacement(pattern, newexpr, idx)

    def sel(self, expr):
        return self._doreplacement(self.varregex, expr, 0)

    def wt(self, expr):
        return self._doreplacement(self.varregex, expr, 1)


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


class CodeLines:
    def __init__(self):
        self.res = []
        self.idlevel = 0

    def Add(self, line):
        self.res.append('  ' * self.idlevel + line)


def GenTreeCode(tree, objlist, multithread=4):
    fname = 'RunTree%i' % GenTreeCode.counter
    # static variable to keep ensure we can give
    # each function a unique name
    GenTreeCode.counter += 1

    # Figure out the list of variable names
    varnames = set()
    varregex = re.compile('[_a-zA-Z][_a-zA-Z0-9]*')
    for obj in objlist:
        varnames.update(varregex.findall(obj.sel))
        for var in obj.var:
            varnames.update(varregex.findall(var))
        varnames.update(varregex.findall(obj.wt))

    code = CodeLines()
    if multithread:
        code.Add('R__LOAD_LIBRARY(libTreePlayer)')
    code.Add('void %s(TTree *tree, TObjArray *hists) {' % fname)
    code.idlevel += 1
    for i, obj in enumerate(objlist):
        code.Add('%s* hist%i = static_cast<%s*>(hists->UncheckedAt(%i));' % (obj.classname, i, obj.classname, i))
    if multithread:
        code.Add('ROOT::EnableImplicitMT(%i);' % multithread)
        for i, obj in enumerate(objlist):
            code.Add('ROOT::TThreadedObject<%s> t_hist%i(*hist%i);' % (obj.classname, i, i))
        code.Add('ROOT::TTreeProcessorMT tp(*tree);')
        code.Add('auto myFunction = [&](TTreeReader &reader) {')
        code.idlevel += 1
    else:
        code.Add('TTreeReader reader(tree);')

    branches = tree.GetListOfBranches()
    binfos = []
    for i in range(branches.GetEntriesFast()):
        branch = branches.UncheckedAt(i)
        bname = branch.GetName()
        btype = branch.GetListOfLeaves().UncheckedAt(0).GetTypeName()
        binfos.append((bname, btype))
        if bname in varnames:
            code.Add('TTreeReaderValue<%s> %s_(reader, "%s");' % (btype, bname, bname))
    if multithread:
        for i, obj in enumerate(objlist):
            code.Add('auto l_hist%i = t_hist%i.Get();' % (i, i))
    code.Add('while (reader.Next()) {')
    code.idlevel += 1
    for bname, btype in binfos:
        if bname in varnames:
            code.Add('%s %s = *%s_;' % (btype, bname, bname))
    for i, obj in enumerate(objlist):
        if multithread:
            code.Add('if (%s) l_hist%i->Fill(%s, %s);' % (obj.sel, i, ', '.join(obj.var), obj.wt))
        else:
            code.Add('if (%s) hist%i->Fill(%s, %s);' % (obj.sel, i, ', '.join(obj.var), obj.wt))

    code.idlevel -= 1
    code.Add('}')
    if multithread:
        code.idlevel -= 1
        code.Add('};')
        code.Add('tp.Process(myFunction);')
        for i, obj in enumerate(objlist):
            code.Add('t_hist%i.Merge()->Copy(*hist%i);' % (i, i))
        code.Add('ROOT::DisableImplicitMT();')
    code.idlevel -= 1
    code.Add('}')
    fullcode = '\n'.join(code.res)
    # print fullcode
    start = time.time()
    ROOT.gInterpreter.Declare(fullcode)
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


def NormaliseTo(hist, val=1.0):
    hist.Scale(val / hist.Integral())
