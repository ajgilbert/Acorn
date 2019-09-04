from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import ROOT


class WGModel(PhysicsModel):
    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("ptBins="):
                self.ptBins = int(po.replace("ptBins=", ""))
            if po.startswith("phiBins="):
                self.phiBins = int(po.replace("phiBins=", ""))
            if po.startswith("type="):
                self.type = str(po.replace("type=", ""))
            if po.startswith("files="):
                self.files = [tuple(X.split(':')) for X in po.replace("files=", "").split(',')]

    def doConstVar(self, name, spec):
        self.modelBuilder.doVar('%s%s' % (name, spec))
        self.modelBuilder.out.var('%s' % name).setConstant(True)

    def doParametersOfInterest(self):
        self.years = ['2016', '2017', '2018']
        self.channels = ['e', 'm']
        self.signs = ['x']
        """Create POI and other parameters, and define the POI set."""
        if self.type in ['eft']:
            self.signs = ['p', 'n']
            # Could measure c3w separately in:
            # fully combined
            # 2016, 2017, 2018
            # electron channel, muon channel
            # W+, W-
            self.doConstVar("c3w", "[0,-10,10]")
            for yr in self.years:
                self.doConstVar("c3w_%s" % yr, "[0,-10,10]")
            for chn in self.channels:
                self.doConstVar("c3w_%s" % chn, "[0,-10,10]")
            for sgn in self.signs:
                self.doConstVar("c3w_%s" % sgn, "[0,-10,10]")
            pois = ['c3w']
            for filename, label, sgn in self.files:
                print filename, label, sgn
                self.ImportFile(filename, sign=sgn, sublabel=label)
            # self.modelBuilder.out.Print()
        elif self.type in ['pt_diff']:
            pois = []
            for sgn in self.signs:
                for ptbin in range(self.ptBins):
                        varname = "r_%s_%i" % (sgn, ptbin)
                        pois.append(varname)
                        self.modelBuilder.doVar("%s[1,0,10]" % varname)
        elif self.type in ['pt_phi_diff']:
            pois = []
            for sgn in self.signs:
                for ptbin in range(self.ptBins):
                    for phibin in range(self.phiBins):
                        varname = "r_%s_%i_%i" % (sgn, ptbin, phibin)
                        pois.append(varname)
                        self.modelBuilder.doVar("%s[1,0,10]" % varname)
        self.modelBuilder.doSet("POI", ",".join(pois))

    def ImportFile(self, filename, sign, sublabel=''):
        fin = ROOT.TFile(filename)
        w = fin.Get('w')
        # w.Print()
        for yr in self.years:
            for chn in self.channels:
                for c in [sign]:
                    self.modelBuilder.factory_('sum::c3w_scale_%s_%s_%s(c3w,c3w_%s,c3w_%s,c3w_%s)' % (yr, chn, c, yr, chn, c))
                    for ptbin in range(self.ptBins):
                        for phibin in range(self.phiBins):
                            name = '%s_%s_%i_%i' % (sublabel, c, ptbin, phibin)
                            # gr = fin.Get(name)
                            # print 'scale_%s_%s_%i_%i' % (sublabel, c, ptbin, phibin)
                            obj = w.function(name)
                            # obj.SetName('scale_%s_%s_%i_%i' % (sublabel, c, ptbin, phibin))
                            # spline = ROOT.RooSpline1D('scale_%s%s_%i_%i' % (c, sub, ptbin, phibin), '', self.modelBuilder.out.var('c3w'), gr.GetN(), gr.GetX(), gr.GetY(), 'CSPLINE')
                            self.modelBuilder.out._import(obj,
                                                          ROOT.RooFit.RenameVariable(name, 'scale_%s_%s_%s_%s_%i_%i' % (yr, sublabel, chn, c, ptbin, phibin)),
                                                          ROOT.RooFit.RenameVariable('c3w', 'c3w_scale_%s_%s_%s' % (yr, chn, c)),
                                                          ROOT.RooFit. RecycleConflictNodes())
        fin.Close()

    def getYieldScale(self, bin, process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"
        if self.DC.isSignal[process]:
            if self.type in ['eft']:
                sgn, chn, ptbin, yr = bin.split('_')
                proc, region, sgn, ptbin, phibin = process.split('_')
                scaler = 'scale_%s_%s_%s_%s_%s_%s' % (yr, region, chn, sgn, ptbin, phibin)
            elif self.type in ['pt_diff']:
                sgn, chn, ptbin, yr = bin.split('_')
                p_args = process.split('_')
                if len(p_args) == 4:
                    proc, region, sgn, ptbin = p_args
                if len(p_args) == 5:
                    proc, region, sgn, ptbin, phibin = p_args
                scaler = 'r_%s_%s' % (sgn, ptbin)
            elif self.type in ['pt_phi_diff']:
                sgn, chn, ptbin, yr = bin.split('_')
                proc, region, sgn, ptbin, phibin = process.split('_')
                scaler = 'r_%s_%s_%s' % (sgn, ptbin, phibin)
            print 'Scaling %s/%s by %s' % (bin, process, scaler)
            return scaler
        return 1

    def getChannelMask(self, bin):
        sgn, chn, ptbin, yr = bin.split('_')
        name = 'mask_%s' % ptbin
        if not self.modelBuilder.out.arg(name):
            self.modelBuilder.doVar('%s[0]' % name)
        return name


wgModel = WGModel()
