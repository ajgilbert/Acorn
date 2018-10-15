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
        self.varregex = re.compile('\\$[_a-zA-Z][_a-zA-Z0-9]*')

    def __getitem__(self, key):
        if key in self.d:
            return self.d[key]
        else:
            self.d[key] = Node()
            return self.d[key]

    def __setitem__(self, key, val):
        self.storage[key] = val

    def _getreplacement(self, matchgroup, override):
        label = matchgroup.group(0)[1:]
        if label in override:
            result = override[label]
        else:
            result = self.storage[label]
        result = '(' + result + ')'
        return result

    def _printMarkers(self, pattern, expr, override):
        verbpositions = []
        substrings = []
        verb_line = ''
        for match in pattern.finditer(expr):
            substrings.append(self._getreplacement(match, override))
            verb_line += ' ' * (match.start(0) - len(verb_line))
            verbpositions.append(len(verb_line))
            verb_line += '^' + '~' * (match.end(0) - match.start(0) - 1)
        if verb_line != '':
            print '\033[1;32;40m' + verb_line + '\033[m'
        for pos, substr in zip(verbpositions, substrings):
            newline = (' ' * pos) + '|' + substr
            printed_newline = (' ' * pos) + '\033[1;32;40m|\033[m'+ substr
            print printed_newline
            self._printMarkers(pattern, newline, override)

    def _doreplacement(self, pattern, expr, override):
        newexpr = pattern.sub(lambda x: self._getreplacement(x, override), expr)
        if newexpr == expr:
            return newexpr
        else:
            return self._doreplacement(pattern, newexpr, override)

    def get(self, expr, override={}, printlevel=0):
        res = self._doreplacement(self.varregex, expr, override)
        if printlevel == 1:
            print '%s --> %s' % (expr, res)
        if printlevel == 2:
            print '%s --> %s' % (expr, res)
            print expr
            self._printMarkers(self.varregex, expr, override)
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
    """
    Args:
        tree (TTree): TTree to process
        objlist (list): List of histograms (TH1, TH2 etc) to fill
        multithread (int): Number of threads to use (0 = no multithreading)

    Returns:
        str: Name of the function that was generated
    """
    fname = 'RunTree%i' % GenTreeCode.counter
    # static variable to keep ensure we can give
    # each function a unique name
    GenTreeCode.counter += 1

    objs_by_sel = defaultdict(list)
    indexed_wts = OrderedDict()

    # Figure out the list of variable names
    varnames = set()
    varregex = re.compile('[_a-zA-Z][_a-zA-Z0-9]*')
    for i, obj in enumerate(objlist):
        obj._draw_idx = i
        varnames.update(varregex.findall(obj.sel))
        for var in obj.var:
            varnames.update(varregex.findall(var))
        varnames.update(varregex.findall(obj.wt))
        objs_by_sel[obj.sel].append(obj)
        # Give each object an index for the weight expression
        if obj.wt not in indexed_wts:
            indexed_wts[obj.wt] = len(indexed_wts)
        obj._wt_idx = indexed_wts[obj.wt]

    code = CodeLines()
    if multithread:
        code.Add('R__LOAD_LIBRARY(libTreePlayer)')
    code.Add('void %s(TTree *tree, TObjArray *hists) {' % fname)
    code.idlevel += 1
    for obj in objlist:
        code.Add('%s* hist%i = static_cast<%s*>(hists->UncheckedAt(%i));' % (obj.classname, obj._draw_idx, obj.classname, obj._draw_idx))
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
            code.Add('std::shared_ptr<%s> l_hist%i = t_hist%i.Get();' % (obj.classname, i, i))
    code.Add('while (reader.Next()) {')
    code.idlevel += 1
    for bname, btype in binfos:
        if bname in varnames:
            code.Add('%s %s = *%s_;' % (btype, bname, bname))
    for wt_expr, wt_idx in indexed_wts.iteritems():
        code.Add('double weightexpr_%i_ = %s;' % (wt_idx, wt_expr))
    for sel, objs in objs_by_sel.iteritems():
        code.Add('if (%s) {' % sel)
        code.idlevel += 1
        for obj in objs:
            if multithread:
                code.Add('l_hist%i->Fill(%s, weightexpr_%i_);' % (obj._draw_idx, ', '.join(obj.var), obj._wt_idx))
            else:
                code.Add('hist%i->Fill(%s, weightexpr_%i_);' % (obj._draw_idx, ', '.join(obj.var), obj._wt_idx))
        code.idlevel -= 1
        code.Add('}')
    code.idlevel -= 1
    code.Add('}')
    if multithread:
        code.idlevel -= 1
        code.Add('};')
        code.Add('tp.Process(myFunction);')
        for obj in objlist:
            code.Add('t_hist%i.Merge()->Copy(*hist%i);' % (obj._draw_idx, obj._draw_idx))
        code.Add('ROOT::DisableImplicitMT();')
    code.idlevel -= 1
    code.Add('}')
    fullcode = '\n'.join(code.res)
    # print fullcode
    start = time.time()
    ROOT.gInterpreter.Declare(fullcode)
    end = time.time()
    print '>> JIT compiled function %s with %i objects in %.2g seconds (%i cores)' % (fname, len(objlist), (end - start), multithread)
    return fname


GenTreeCode.counter = 0


def MultiDraw(node, sample_to_file_dict, tree_name, mt_cores=0, mt_thresh=1E7):
    drawtree = defaultdict(list)
    codegen_time = 0
    draw_time = 0
    for path, name, obj in node.ListObjects():
        drawtree[obj.sample].append(obj)
    # pprint(drawtree)
    for sample, drawobjs in drawtree.iteritems():
        f = ROOT.TFile(sample_to_file_dict[sample])
        t = f.Get(tree_name)
        entries = t.GetEntriesFast()
        use_cores = 0
        if mt_cores > 0 and entries > mt_thresh:
            use_cores = mt_cores
        start = time.time()
        fname = GenTreeCode(t, drawobjs, use_cores)
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
    if hist.Integral() != 0.0:
        hist.Scale(val / hist.Integral())

def WidthDivide(hist):
        hist.Scale(1., 'width')
