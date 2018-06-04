import ROOT
import sys

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)

if len(sys.argv) < 3:
    sys.exit(1)

fin = ROOT.TFile(sys.argv[1].split(':')[0])

lumiobj = fin.Get(sys.argv[1].split(':')[1])

with open(sys.argv[2], "w") as text_file:
    text_file.write(lumiobj.AsJsonString())
