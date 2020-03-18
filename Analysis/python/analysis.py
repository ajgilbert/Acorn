import ROOT
import time
import re
import math
import hashlib
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

    def ListNodes(self, depth=9999999, res=None, dirname='', withObjects=False):
        if res is None:
            res = []
        pre = '/' if dirname is not '' else ''
        if withObjects:
            hasObjects = False
            for key, val in self.d.iteritems():
                if not isinstance(val, Node):
                    hasObjects = True
                    break
            if hasObjects:
                res.append((dirname, self))
        else:
            res.append((dirname, self))
        if depth > 0:
            for key, val in self.d.iteritems():
                if isinstance(val, Node):
                    val.ListNodes(depth - 1, res, dirname + pre + key, withObjects)
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


def TDirToNode(file, startdir='', node=None, onlyDirs=False):
    if node is None:
        node = Node()
    file.cd()
    ROOT.gDirectory.cd(startdir + '/')
    dirs_to_proc = []
    for key in ROOT.gDirectory.GetListOfKeys():
        if key.GetClassName() == 'TDirectoryFile':
            dirs_to_proc.append(key.GetName())
        elif not onlyDirs:
            obj = key.ReadObj()
            node[obj.GetName()] = obj
    for dirname in dirs_to_proc:
        TDirToNode(file, startdir + '/' + dirname, node[dirname], onlyDirs)
    return node


def NodeToTDir(file, node, postfix=None):
    post = ''
    if postfix is not None:
        post = postfix
    startdir = ROOT.gDirectory.GetPath()
    for path, subnode in node.ListNodes():
        file.cd()
        as_vec = path.split('/')
        if len(as_vec) >= 1:
            for i in xrange(0, len(as_vec)):
                if not ROOT.gDirectory.GetDirectory(as_vec[i]):
                    ROOT.gDirectory.mkdir(as_vec[i])
                ROOT.gDirectory.cd(as_vec[i])
        for opath, name, obj in subnode.ListObjects(depth=0):
            if not ROOT.gDirectory.FindKey(name + post):
                obj.SetName(name + post)
                ROOT.gDirectory.WriteTObject(obj, name + post)
        ROOT.gDirectory.cd(startdir)


def GetListOfDirectories(file, dirlist=None, curr=[]):
    if dirlist is None:
        dirlist = list()
        file.cd()
    ROOT.gDirectory.cd('/' + '/'.join(curr))
    dirs_to_proc = []
    has_other_objects = False
    for key in ROOT.gDirectory.GetListOfKeys():
        if key.GetClassName() == 'TDirectoryFile':
            dirs_to_proc.append(key.GetName())
        else:
            has_other_objects = True
    if has_other_objects:
        # This directory should be added to the list
        dirlist.append('/' + '/'.join(curr))
    for dirname in dirs_to_proc:
        GetListOfDirectories(file, dirlist, curr + [dirname])
    return dirlist


def GetHistsInDir(file, dir):
    res = {}
    file.cd()
    ROOT.gDirectory.cd(dir)
    for key in ROOT.gDirectory.GetListOfKeys():
        if key.GetClassName() != 'TDirectoryFile':
            obj = key.ReadObj()
            res[obj.GetName()] = obj
    return res

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
    # res.SetTitle(json.dumps({
    #     'sample': sample,
    #     'var': var,
    #     'wt': wt,
    #     'sel': sel
    #     }))
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

    classnames = list(objlist.keys())
    nobjs_total = 0
    for clname in classnames:
        nobjs_total += len(objlist[clname])
        for i, obj in enumerate(objlist[clname]):
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
    obj_offset = 0  # offset in the full list
    for clname in classnames:
        nobjs = len(objlist[clname])
        code.Add('std::vector<%s *> v_%s_hists(%i);' % (clname, clname, nobjs))
        code.Add('for (unsigned i = 0; i < %i; ++i) {' % nobjs)
        code.idlevel += 1
        code.Add('v_%s_hists[i] = static_cast<%s*>(hists->UncheckedAt(i + %i));' % (clname, clname, obj_offset))
        code.idlevel -= 1
        code.Add('}')
        obj_offset += nobjs
    if multithread:
        code.Add('ROOT::EnableImplicitMT(%i);' % multithread)
        for clname in classnames:
            nobjs = len(objlist[clname])
            code.Add('std::vector<std::unique_ptr<ROOT::TThreadedObject<%s>>> v_%s_t_hists(%i);' % (clname, clname, nobjs))
            code.Add('for (unsigned i = 0; i < %i; ++i) {' % nobjs)
            code.idlevel += 1
            code.Add('v_%s_t_hists[i] = std::make_unique<ROOT::TThreadedObject<%s>>(*(v_%s_hists[i]));' % (clname, clname, clname))
            code.idlevel -= 1
            code.Add('}')
        code.Add('ROOT::TTreeProcessorMT tp(*tree);')
        code.Add('unsigned ncalled = 0;')
        code.Add('auto myFunction = [&](TTreeReader &reader) {')
        code.idlevel += 1
        code.Add('ncalled += 1;')
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
        for clname in classnames:
            nobjs = len(objlist[clname])
            code.Add('std::vector<std::shared_ptr<%s>> v_%s_l_hists(%i);' % (clname, clname, nobjs))
            code.Add('for (unsigned i = 0; i < %i; ++i) {' % nobjs)
            code.idlevel += 1
            code.Add('v_%s_l_hists[i] = v_%s_t_hists[i]->Get();' % (clname, clname))
            code.idlevel -= 1
            code.Add('}')

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
                code.Add('v_%s_l_hists[%i]->Fill(%s, weightexpr_%i_);' % (obj.classname, obj._draw_idx, ', '.join(obj.var), obj._wt_idx))
            else:
                code.Add('v_%s_hists[%i]->Fill(%s, weightexpr_%i_);' % (obj.classname, obj._draw_idx, ', '.join(obj.var), obj._wt_idx))
        code.idlevel -= 1
        code.Add('}')
    code.idlevel -= 1
    code.Add('}')
    if multithread:
        code.idlevel -= 1
        code.Add('};')
        # Sometimes it seems small Trees don't have proper clusters assigned (we end up splitting with 1 per event)
        # This seems to be the only way to detect this...
        code.Add('auto it = tree->GetClusterIterator(0); it.Next();')
        code.Add('if (it.GetNextEntry() == 1) {')
        code.idlevel += 1
        code.Add('std::cout << ">> Tree does not have proper clusters, running single-threaded" << std::endl;')
        code.Add('TTreeReader reader(tree);')
        code.Add('myFunction(reader);')
        code.idlevel -= 1
        code.Add('} else {')
        code.idlevel += 1
        code.Add('tp.Process(myFunction);')
        code.idlevel -= 1
        code.Add('}')
        for clname in classnames:
            nobjs = len(objlist[clname])
            code.Add('for (unsigned i = 0; i < %i; ++i) {' % nobjs)
            code.idlevel += 1
            code.Add('v_%s_t_hists[i]->Merge()->Copy(*(v_%s_hists[i]));' % (clname, clname))
            code.idlevel -= 1
            code.Add('}')
        code.Add('ROOT::DisableImplicitMT();')
        code.Add('std::cout << ">> Draw function was called " << ncalled << " times" << std::endl;')
    code.idlevel -= 1
    code.Add('}')
    fullcode = '\n'.join(code.res)
    # If we want to dump the code for testing...
    # includes = [
    #     '#include "TTreeReader.h"',
    #     '#include "TTree.h"',
    #     '#include "TObjArray.h"',
    #     '#include "TH1F.h"',
    #     '#include "TH2F.h"',
    #     '#include "ROOT/TTreeProcessorMT.hxx"',
    #     '#include <memory>',
    #     ''
    # ]
    # with open('%s.cc' % fname, 'w') as out_file:
    #     out_file.write('\n'.join(includes + code.res))
    # print fullcode
    start = time.time()
    ROOT.gInterpreter.Declare(fullcode)
    end = time.time()
    print '>> JIT compiled function %s with %i objects in %.2g seconds (%i cores)' % (fname, nobjs_total, (end - start), multithread)
    return fname


GenTreeCode.counter = 0


def DrawObjsHash(drawobjs):
    allargs = []
    for obj in drawobjs:
        allargs.extend([obj.classname, obj.var, obj.sel, obj.wt])
    return hashlib.md5(str(allargs)).hexdigest()


def MultiDraw(node, sample_to_file_dict, tree_name, mt_cores=0, mt_thresh=1E7):
    drawtree = defaultdict(lambda: defaultdict(list))
    codegen_time = 0
    draw_time = 0
    for path, name, obj in node.ListObjects():
        drawtree[obj.sample][obj.classname].append(obj)
    # pprint(drawtree)
    # sys.exit(0)
    func_dict = {}
    for sample, drawobjs in drawtree.iteritems():
        flatobjlist = []
        for objlist in drawobjs.values():
            flatobjlist.extend(objlist)

        f = ROOT.TFile.Open(sample_to_file_dict[sample])
        t = f.Get(tree_name)
        entries = t.GetEntriesFast()
        use_cores = 0
        if mt_cores > 0 and entries > mt_thresh:
            use_cores = mt_cores
        start = time.time()
        obj_hash = DrawObjsHash(flatobjlist)
        if obj_hash not in func_dict:
            fname = GenTreeCode(t, drawobjs, use_cores)
            func_dict[obj_hash] = fname
        else:
            fname = func_dict[obj_hash]
            print '>> Recycling function %s in %.2g seconds' % (fname, (time.time() - start))
        end = time.time()
        codegen_time += (end - start)
        histarr = ROOT.TObjArray(len(flatobjlist))
        for i, h in enumerate(flatobjlist):
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


def VariableRebin(hist, binning):
    newhist = hist.Rebin(len(binning) - 1, "", array('d', binning))
    return newhist


def ReadTxtHist(filename, column=1, column_err=2, classname='TH1D', objname='hist', correct_density=True):
    with open(filename) as infile:
        lines = [tuple(line.split()) for line in infile.readlines() if not line.startswith('#')]
    N = (len(lines) - 1) / 2
    edges = []
    contents = []
    errors = []
    for i in xrange(N):
        lo_edge = float(lines[(i * 2)][0])
        hi_edge = float(lines[(i * 2) + 1][0])
        edges.append(lo_edge)
        if i == N - 1:
            edges.append(hi_edge)
        width = 1.0
        if correct_density:
            width = hi_edge - lo_edge
        contents.append(float(lines[i * 2][column]) * width)
        errors.append(float(lines[i * 2][column_err]))
    hist = getattr(ROOT, classname)(objname, objname, len(contents), array('d', edges))
    for i in xrange(N):
        hist.SetBinContent(i + 1, contents[i])
        hist.SetBinError(i + 1, errors[i])
    return hist
    print hist.Integral()
    print hist.Integral('width')


def BinningFromStr(binning_str):
    if binning_str.startswith('['):
        binning = [float(eval(v)) for v in binning_str[1:-1].split(',')]
        binning = (len(binning) - 1, array('d', binning))
    else:
        binning = binning_str[1:-1].split(',')
        binning = (int(binning[0]), float(eval(binning[1])), float(eval(binning[2])))
    return binning


def StrFromBinning(bins):
    if isinstance(bins, list):
        return '[%s]' % ','.join([str(x) for x in bins])
    if isinstance(bins, tuple):
        return '(%i,%f,%f)' % bins


def BinEdgesFromStr(binning_str):
    binning = BinningFromStr(binning_str)
    if binning_str.startswith('['):
        return list(binning[1])
    else:
        width = (float(binning[2]) - float(binning[1])) / float(binning[0])
        binedges = []
        for i in range(binning[0]):
            binedges.append(float(binning[1]) + float(i) * width)
        binedges.append(binning[2])
        return binedges


def SetupPads(split_points, gaps_left, gaps_right, ratio=None, purity=None):
    pads = []
    ratio_pads = []
    purity_pads = []
    l_margin = ROOT.gStyle.GetPadLeftMargin()
    r_margin = ROOT.gStyle.GetPadRightMargin()
    t_margin = ROOT.gStyle.GetPadTopMargin()
    b_margin = ROOT.gStyle.GetPadBottomMargin()
    usable_width = 1. - l_margin - r_margin
    usable_height = 1. - t_margin - b_margin
    for i in xrange(len(split_points)+1):
        pad = ROOT.TPad('pad%i'%i, '', 0., 0., 1., 1.)
        if i > 0:
            pad_l_margin = l_margin + sum(split_points[0:i]) * usable_width + gaps_left[i-1]
            pad.SetLeftMargin(pad_l_margin)
        if i < len(split_points):
            pad_r_margin = (1. - sum(split_points[0:i+1])) * usable_width + r_margin + gaps_right[i]
            pad.SetRightMargin(pad_r_margin)
        print pad.GetLeftMargin(), pad.GetRightMargin()
        if ratio is not None:
            if purity is None:
                pad.SetBottomMargin(b_margin + usable_height * ratio)
            else:
                pad.SetBottomMargin(b_margin + usable_height * (ratio + purity))

        pad.SetFillStyle(4000)
        pad.Draw()
        pads.append(pad)
        if ratio is not None:
            padr = ROOT.TPad('r_pad%i'%i, '', 0., 0., 1., 1.)
            if i > 0:
                padr.SetLeftMargin(pad_l_margin)
            if i < len(split_points):
                padr.SetRightMargin(pad_r_margin)
            padr.SetTopMargin(1 - (b_margin + usable_height * ratio))
            padr.SetFillStyle(4000)
            padr.Draw()
            ratio_pads.append(padr)
        if purity is not None:
            padp = ROOT.TPad('p_pad%i'%i, '', 0., 0., 1., 1.)
            if i > 0:
                padp.SetLeftMargin(pad_l_margin)
            if i < len(split_points):
                padp.SetRightMargin(pad_r_margin)
            padp.SetTopMargin(1 - (b_margin + usable_height * (ratio + purity)))
            padp.SetBottomMargin(b_margin + usable_height * ratio)
            padp.SetFillStyle(4000)
            padp.Draw()
            purity_pads.append(padp)
    pads[0].cd()
    # pads.reverse()
    return pads, ratio_pads, purity_pads


def CheckBinErrors(node, hists):
    N = node[hists[0]].GetNbinsX()
    for i in xrange(N):
        toterr = 0.0
        totval = 0.0
        errs = []
        for h in hists:
            err = node[h].GetBinError(i + 1)
            val = node[h].GetBinContent(i + 1)
            errs.append((h, err, val))
            totval += val
            toterr += (err * err)
        print '>> Bin %i has total error %f / %f' % (i, math.sqrt(toterr), totval)
        print sorted(errs, key=lambda x: x[1], reverse=True)[:5]
