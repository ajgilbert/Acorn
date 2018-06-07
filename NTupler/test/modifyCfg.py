import argparse
import os
import json
import sys
import math

edits = []

class Visitor:
    def __init__(self, modtype):
        self.modtype = modtype
        self.results = []

    def enter(self, obj):
        if isinstance(obj, cms.OutputModule):
            self.results.append(obj)

    def leave(self, obj):
        pass

def UpdateModuleIfExists(mod, attrName, newVal, prefix='process.'):
    if hasattr(mod, attrName):
        attr = getattr(mod, attrName)
        prev = '%s' % attr
        attr.setValue(newVal)
        print 'Updating %s from %s to %s' % (attrName, prev, getattr(mod, attrName))
        edits.append('%s%s.%s = %s' % (prefix, mod, attrName, getattr(mod, attrName)))

def UpdateIfExists(process, modName, attrName, newVal, prefix='process.'):
    if hasattr(process, modName):
        mod = getattr(process, modName)
        if hasattr(mod, attrName):
            attr = getattr(mod, attrName)
            prev = '%s' % attr
            newArg = type(attr)(newVal)
            attr.setValue(newVal)
            print 'Updating %s.%s from %s to %s' % (modName, attrName, prev, getattr(mod, attrName))
            edits.append('%s%s.%s = %s' % (prefix, modName, attrName, getattr(mod, attrName)))

def Update(process, modName, attrName, newVal, prefix='process.'):
    if hasattr(process, modName):
        mod = getattr(process, modName)
        setattr(mod, attrName, newVal)
        print 'Updating %s.%s to %s' % (modName, attrName, getattr(mod, attrName))
        edits.append('%s%s.%s = %s' % (prefix, modName, attrName, getattr(mod, attrName)))

def SetRandomSeeds(process, newSeed):
    if not hasattr(process, 'RandomNumberGeneratorService'):
        return
    mod = getattr(process, 'RandomNumberGeneratorService')
    for par in mod.parameterNames_():
        UpdateIfExists(mod, par, 'initialSeed', newSeed, prefix='process.RandomNumberGeneratorService.')

def SetEvents(process, nEvents):
    UpdateIfExists(process, 'externalLHEProducer', 'nEvents', nEvents)
    UpdateIfExists(process, 'maxEvents', 'input', nEvents)

def FindOutputModule(process):
    visitor = Visitor(cms.OutputModule)
    for schd in process.schedule:
        schd.visit(visitor)
    if len(visitor.results) == 1:
        return visitor.results[0]
    else:
        raise RuntimeError('Too many OutputModules!')

def SetOutputFileName(process, newName, namedModule=None):
    if namedModule is not None:
        UpdateModuleIfExists(getattr(process, namedModule), 'fileName', 'file:%s' % newName)
    else:
        UpdateModuleIfExists(FindOutputModule(process), 'fileName', 'file:%s' % newName)
def SetInputFileName(process, newName):
    UpdateIfExists(process, 'source', 'fileNames', ['file:%s' % newName])


def UpdateConfig(inputCfg, outputCfg, events=None, randomSeeds=None, inputFile=None, outputFile=None, outputModule=None, setLumiOffsets=None):
    # Have to reset sys.argv here in case the inputCfg will do its own VarParsing
    # => this is a bit of a hack!
    sys.argv = [inputCfg]
    execfile(inputCfg, globals())
    if events is not None:
        SetEvents(process, int(events))
    if randomSeeds is not None:    
        SetRandomSeeds(process, int(randomSeeds))
    if inputFile is not None:
        SetInputFileName(process, str(inputFile))
    if outputFile is not None:
        SetOutputFileName(process, str(outputFile), namedModule=outputModule)

    if setLumiOffsets is not None: 
        eventsPerLumi = int(setLumiOffsets)
        firstLumi = 1 + (int(randomSeeds) - 1) * int(math.ceil(float(events) / float(eventsPerLumi)))
        Update(process, 'source', 'numberEventsInLuminosityBlock', cms.untracked.uint32(eventsPerLumi))
        Update(process, 'source', 'firstLuminosityBlock', cms.untracked.uint32(firstLumi))

    if args.strategy == 0:
        outFile = open(outputCfg, "w")
        outFile.write(process.dumpPython())
        outFile.close()
    elif args.strategy == 1:
        os.system('cp %s %s' % (inputCfg, outputCfg))
        with open(outputCfg) as ofile:
            cfg = ofile.read()
        cfg = cfg.replace('#{EDITS}', '\n'.join(edits))
        outFile = open(outputCfg, "w")
        outFile.write(cfg)
        outFile.close()



parser = argparse.ArgumentParser()
parser.add_argument('io', nargs=2,
                    help='[input cfg] [output cfg]')
parser.add_argument('--events', type=int, default=None, help='Number of events to process')
parser.add_argument('--randomSeeds', type=int, default=None, help='Set random seeds')
parser.add_argument('--inputFile', type=str, default=None, help='Set input file')
parser.add_argument('--outputFile', type=str, default=None, help='Set output file')
parser.add_argument('--strategy', type=int, default=1, help='Patching strategy')
parser.add_argument('--setLumiOffsets', type=int, default=None, help='Set this many events per lumiBlock')
parser.add_argument('--outputModule', type=str, default=None, help='Output module to modify')

args = parser.parse_args()


UpdateConfig(args.io[0], args.io[1], events=args.events, randomSeeds=args.randomSeeds, inputFile=args.inputFile, outputFile=args.outputFile, outputModule=args.outputModule, setLumiOffsets=args.setLumiOffsets)

