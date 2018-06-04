import os

paths = [
    'HLT_IsoMu24_v',
    'HLT_IsoTkMu24_v',
    'HLT_IsoMu27_v'
    ]

for path in paths:
    print 'brilcalc trg --prescale --hltpath "%s*" -o %s.csv' % (path, path)
    os.system('brilcalc trg --prescale --hltpath "%s*" -o %s.csv' % (path, path))
