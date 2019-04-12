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

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.type in ['eft']:
            self.modelBuilder.doVar("c3w[0,0,10]")
            pois = ['c3w']
            for filename, label in self.files:
                print filename, label
                self.ImportFile(filename, 'p', sublabel=label)
        if self.type in ['pt_diff']:
            pois = []
            for c in ['p', 'n']:
                for ptbin in range(self.ptBins):
                        varname = "r_%s_%i" % (c, ptbin)
                        pois.append(varname)
                        self.modelBuilder.doVar("%s[1,0,10]" % varname)
        else:
            pois = []
            for c in ['p', 'n']:
                for ptbin in range(4):
                    for phibin in range(5):
                        varname = "r_%s_%i_%i" % (c, ptbin, phibin)
                        pois.append(varname)
                        self.modelBuilder.doVar("%s[1,-10,10]" % varname)
        self.modelBuilder.doSet("POI", ",".join(pois))

    def ImportFile(self, filename, sign, sublabel=''):
        fin = ROOT.TFile(filename)
        w = fin.Get('w')
        w.Print()
        for c in [sign]:
            for ptbin in range(self.ptBins):
                for phibin in range(self.phiBins):
                    name = '%s_w_%s_bin_%i_%i' % (sublabel, c, ptbin + 1, phibin + 1)
                    # gr = fin.Get(name)
                    sub = ''
                    if sublabel != 'main':
                        sub = str('_' + sublabel)
                    print 'scale_%s%s_%i_%i' % (c, sub, ptbin, phibin)
                    obj = w.function(name)
                    obj.SetName('scale_%s%s_%i_%i' % (c, sub, ptbin, phibin))
                    # spline = ROOT.RooSpline1D('scale_%s%s_%i_%i' % (c, sub, ptbin, phibin), '', self.modelBuilder.out.var('c3w'), gr.GetN(), gr.GetX(), gr.GetY(), 'CSPLINE')
                    self.modelBuilder.out._import(obj)
        fin.Close()

    def getYieldScale(self, bin, process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"
        if self.DC.isSignal[process]:
            if self.type in ['eft']:
                scaler = process.replace('WG_', 'scale_')
            elif self.type in ['pt_diff']:
                scaler = process.replace('WG_', 'r_')
                scaler=  scaler.replace('_met1', '')
            print 'Scaling %s/%s by %s' % (bin, process, scaler)
            return scaler
        return 1


wgModel = WGModel()
