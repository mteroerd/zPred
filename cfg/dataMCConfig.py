import sys
import pickle
from frameworkStructure import pathes
from locations import locations
from defs import getRegion, getPlot, getRunRange, Backgrounds

import ROOT
from helpers import readTrees, getDataHist, TheStack, totalNumberOfGeneratedEvents, Process, applyCorrections

cutString =     {       "pos" : "(JZBCONDITION > 0)",
                        "neg" : "(JZBCONDITION < 0)",
                        "2j"  : "(nJets == 2)",
                        "3j"  : "(nJets >= 3)",
                        "zPeak":"(p4.M() < 101 && p4.M() > 81)"
                }

def configurePlot(plot, jzbType, plotData, direction, responseCorr, puCorr, peakCorr,correctMET,correctionMode):
        from defs import varToUse
        fileName = "%s_%d.pkl"%(jzbType,correctionMode)
        with open("shelves/"+fileName, "r") as corrFile:
                corrs = pickle.load(corrFile)
                corrFile.close()
                corrs = corrs[plotData][direction]
                
                for key in varToUse.keys():
                        plot.variable = plot.variable.replace(key, varToUse[key][jzbType])
                        if key == "MET" and (not correctMET):
                                applyCorrections(plot, corrs, False, False, False, mode="var")
                        else:
                                applyCorrections(plot, corrs, responseCorr, puCorr, peakCorr, mode="var")
                plot.cuts = plot.cuts.replace("JZBCONDITION", varToUse["JZB"][jzbType])
                applyCorrections(plot, corrs, responseCorr, puCorr, peakCorr,mode="cut")
                plot.cuts = plot.cuts.replace("METCONDITION", varToUse["MET"][jzbType])
                if correctMET:
                        applyCorrections(plot, corrs, responseCorr, puCorr, peakCorr,mode="cut")
                else:
                        applyCorrections(plot, corrs, False, False, False,mode="cut")
class dataMCConfig:
        #jzbType = "unCorrMet"
        jzbType = "type-IMet"
        
        onlyShift = True
        correctionMode = 1 # 3=response+pu+peak, 2=pu+peak, 1=peak, 0=none
        
        responseCorr = False
        puCorr = False
        peakCorr = False
        correctMET = True
        
        plotData = True
        plotSyst = False
        normalizeToData = False
        plotRatio = True
        plotSignal = False
        useTriggerEmulation = False 
        doPUWeights = False     
        DontScaleTrig = False 
        personalWork = True
        doTopReweighting = True
        preliminary = True
        forPAS = False
        forTWIKI = False
        
        plot = None
        plot2 = None
        region = None
        runRange = None
        dataSetPath = ""
        prefix = ""
        theoUncert = 0.
        backgrounds = []
                
        def __init__(self,plot,plot2=None,region="Inclusive",runName = "Run2015_25ns",plotData=True,normalizeToData=False,plotRatio=True,signals=None,useTriggerEmulation=False,personalWork=False,doTopReweighting=False,preliminary=True,forPAS=False,forTWIKI=False,backgrounds = [],produceTheoUncert=False,dontScaleTrig=False,plotSyst=False,doPUWeights=False,jzbType = "None",responseCorr = False ,puCorr = False, peakCorr = False, correctMET = True):
                sys.path.append(pathes.basePath)
                if jzbType == "None":
                        jzbType = dataMCConfig.jzbType
                
                self.jzbType = jzbType
                self.prefix = jzbType 
                
                self.responseCorr = responseCorr
                self.puCorr = puCorr
                self.peakCorr = peakCorr
                self.correctMET = correctMET
                 
                self.dataSetPath = locations.dataSetPath
                if dontScaleTrig:
                        self.dataSetPath = locations.dataSetPathTrigger
                self.runRange = getRunRange(runName)
                
                self.selection = getRegion(region)
                
                self.plotData = plotData
                
                plotList = None
                if "#" in plot:
                        plotList = plot.split("#")
                        plot = plotList[0]
                        plotList = plotList[1:]
                
                self.plot = getPlot(plot)
                self.plot.addRegion(self.selection)
                self.plot.cleanCuts()
                self.plot.cuts = self.plot.cuts % self.runRange.runCut
                
                if plotList != None:
                        self.plot.cuts = self.plot.cuts[:-1]
                        for l in plotList:
                                self.plot.cuts = self.plot.cuts + "&& %s"%(cutString[l])
                        self.plot.cuts = self.plot.cuts + ")"
                
                
                if plot2 != None:
                        plotList = None
                        if "#" in plot2:
                                plotList = plot2.split("#")
                                plot2 = plotList[0]
                                plotList = plotList[1:]
                        
                        self.plot2 = getPlot(plot2)
                        self.plot2.addRegion(self.selection)
                        self.plot2.cleanCuts()
                        self.plot2.cuts = self.plot2.cuts % self.runRange.runCut
                        
                        if plotList != None:
                                self.plot2.cuts = self.plot2.cuts[:-1]
                                for l in plotList:
                                        self.plot2.cuts = self.plot2.cuts + " && %s"%(cutString[l])
                                self.plot2.cuts = self.plot2.cuts + ")"
                if "Inclusive" in region:
                        self.direction = "Inclusive"
                elif "Central" in region:
                        self.direction = "Central"
                elif "Forward" in region:
                        self.direction = "Forward"
                
                
                configurePlot(self.plot, self.jzbType, self.plotData, self.direction, self.responseCorr, self.puCorr, self.peakCorr, self.correctMET, self.correctionMode)
                if plot2 != None:
                        configurePlot(self.plot2, self.jzbType, self.plotData, self.direction, self.responseCorr, self.puCorr, self.peakCorr, self.correctMET, self.correctionMode)
                
                self.normalizeToData = normalizeToData
                self.plotRatio = plotRatio
                self.signals = signals
                if self.signals is not None:
                        self.plotSignal = True
                self.backgrounds = backgrounds
                self.useTriggerEmulation = useTriggerEmulation
                self.personalWork = personalWork
                self.preliminary = preliminary
                self.doTopReweighting = doTopReweighting
                self.forPAS = forPAS
                self.forTWIKI = forTWIKI
                self.DontScaleTrig = dontScaleTrig
                self.doPUWeights = doPUWeights
                
                from corrections import rSFOF   
                self.rSFOF = rSFOF
                
               
