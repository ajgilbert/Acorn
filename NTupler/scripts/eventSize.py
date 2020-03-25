import ROOT
import os
import json
import argparse
import subprocess
import imp
import re
import pprint

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('-t', '--tree', default='EventTree')
parser.add_argument('-s', '--sub-branches', action='store_true')
args = parser.parse_args()


dosub = args.sub_branches

fin = ROOT.TFile(args.input)
tree = fin.Get(args.tree)
entries = tree.GetEntriesFast()
totbytes = float(tree.GetTotBytes())
zipbytes = float(tree.GetZipBytes())

branches = tree.GetListOfBranches()
results = []
results_sub = []
for i in range(branches.GetSize()):
    branch = branches[i]
    if not branch:
        continue
    res = {
        'Name': branch.GetName(),
        'TotBytes': float(branch.GetTotBytes("*")),
        'ZipBytes': float(branch.GetZipBytes("*")),
    }
    res['Compression'] = res['TotBytes'] / res['ZipBytes']
    res['TotFrac'] = res['TotBytes'] / totbytes
    res['ZipFrac'] = res['ZipBytes'] / zipbytes
    res['ZipKBperEvent'] = (res['ZipBytes'] / 1024) / entries
    results.append(res)

    subbranches = branch.GetListOfBranches()
    for j in range(subbranches.GetSize()):
        subbranch = subbranches[j]
        if not subbranch:
            continue
        res = {
            'Name': subbranch.GetName(),
            'TotBytes': float(subbranch.GetTotBytes("*")),
            'ZipBytes': float(subbranch.GetZipBytes("*")),
        }
        res['Compression'] = res['TotBytes'] / res['ZipBytes']
        res['TotFrac'] = res['TotBytes'] / totbytes
        res['ZipFrac'] = res['ZipBytes'] / zipbytes
        res['ZipKBperEvent'] = (res['ZipBytes'] / 1024) / entries
        results_sub.append(res)


print 'Analysing: %s' % (args.tree)
print '%-20s : %i' % ('Entries', entries)
print '%-20s : %.2f' % ('TotalMB', totbytes / 1048576.0)
print '%-20s : %.2f' % ('ZippedMB', zipbytes / 1048576.0)
print '%-20s : %.2f' % ('Compression', totbytes / zipbytes)
print '%-20s : %.2f' % ('ZipKBperEvent', (zipbytes / 1024.) / entries)
print '%-20s : %.2f' % ('1MilFileSizeMB', (zipbytes / 1048576.0) * (1E6 / float(entries)))

sorted_results = sorted(results, key=lambda x: x['ZipFrac'], reverse=True)
sorted_results_sub = sorted(results_sub, key=lambda x: x['ZipFrac'], reverse=True)

for b in sorted_results:
    print '%-50s : %.1f : %.3f' % (b['Name'], b['ZipFrac'] * 100., b['ZipKBperEvent'])
if dosub:
    print '-------'
    for b in sorted_results_sub:
        print '%-80s : %.1f : %.3f' % (b['Name'], b['ZipFrac'] * 100., b['ZipKBperEvent'])

# print 'Entries: %i, TotalMB: %.2f, ZippedMB: %.2f, Compression factor: %.2f' % (entries, totbytes / 1048576.0, zipbytes / 1048576.0, zipbytes / totbytes)
