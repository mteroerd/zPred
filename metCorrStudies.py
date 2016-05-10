import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np


from ROOT import TMath

import argparse 
import dataMCConfig

import gc
gc.enable() 

from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TMath, gROOT, TF1, TGraphAsymmErrors,gStyle
import ratios
from defs import defineMyColors
from defs import myColors   
from defs import Region
from defs import Regions
from defs import Backgrounds
from defs import Plot, getPlot
from setTDRStyle import setTDRStyle
gROOT.SetBatch(True)
from helpers import *   
import math

from centralConfig import plotLists
from sampleConfig import getBackgrounds
from plotTemplate import plotTemplate2D

def make2DPlot(dilepton, plot1, plot2, plotData=False):
    bkg = getBackgrounds("DY", "TT")
    mainConfig = dataMCConfig.dataMCConfig(plot=plot1,plot2=plot2,region="Inclusive",runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False)
    plott2 = plot2

    eventCounts = totalNumberOfGeneratedEvents(mainConfig.dataSetPath)  

    processes = []
    for background in mainConfig.backgrounds:
        processes.append(Process(getattr(Backgrounds,background),eventCounts))

    template = plotTemplate2D(mainConfig)
    template.labelZ = "Events"
    template.logZ = True
    
    template.dilepton = dilepton
    
    jzbInd = "JZB > 0" if "pos" in plot1 else "JZB < 0"
    jetInd = "#geq3j" if "3j" in plot1 else "=2j"
    template.cutsText = jetInd+", "+jzbInd
    
    treeEE = readTrees(mainConfig.dataSetPath, "EE")
    treeMuMu = readTrees(mainConfig.dataSetPath, "MuMu")
    treeEMu = readTrees(mainConfig.dataSetPath, "EMu")
 
    mainConfig.plot.addDilepton(dilepton)    
    plot = mainConfig.plot
    plot2 = mainConfig.plot2
    
    scaleTree1 = 1.0
    scaleTree2 = 1.0
    if mainConfig.plot.tree1 == "EE":
        tree1 = treeEE
        scaleTree1 = mainConfig.selection.trigEffs.effEE.val
    elif mainConfig.plot.tree1 == "MuMu":
        tree1 = treeMuMu
        scaleTree1 = mainConfig.selection.trigEffs.effMM.val
    elif mainConfig.plot.tree1 == "EMu":
        tree1 = treeEMu 
        scaleTree1 = mainConfig.selection.trigEffs.effEM.val            
    else: 
        print "Unknown Dilepton combination! %s not created!"%(mainConfig.plot.filename,)
        return
    
    if mainConfig.plot.tree2 != "None":
        if mainConfig.plot.tree2 == "EE":
                tree2 = treeEE
                scaleTree2 = mainConfig.selection.trigEffs.effEE.val                
        elif mainConfig.plot.tree2 == "MuMu":
                tree2 = treeMuMu
                scaleTree2 = mainConfig.selection.trigEffs.effMM.val

        elif mainConfig.plot.tree2 == "EMu":
                tree2 = treeEMu 
                scaleTree2 = mainConfig.selection.trigEffs.effEM.val                    
        else:
            print "Unknown Dilepton combination! %s not created!"%(mainConfig.plot.filename,)
            return
    else:
        tree2 = "None"
        
    if mainConfig.useTriggerEmulation or mainConfig.DontScaleTrig:
        scaleTree2 = 1.0
        scaleTree1 = 1.0                
    
    gStyle.SetNumberContours(30)
    if mainConfig.plotData:
        fullHist = getData2DHist(plot,plot2, tree1, tree2)
    else:
        MCHists = []
        for process in processes:
            MCHists.append(process.createCombined2DHistogram(mainConfig.runRange.lumi,plot,plot2,tree1,tree2,1,scaleTree1,scaleTree2,doTopReweighting=mainConfig.doTopReweighting, doPUWeights=mainConfig.doPUWeights))
        fullHist = MCHists[0].Clone()
        fullHist.Reset()
        for hist in MCHists:
            fullHist.Add(hist,1)
    
    template.setPrimaryPlot(fullHist, "COLZ")
    
    template.draw()
    template.useJZBPath = False
    template.setFolderName("metCorrStudies")
    template.saveAs(plot1+"_"+plott2)

    
def main():
    import multiprocessing as mp
    import time
    t1 = time.time()

    processes = []
    #processes += [mp.Process(target=make2DPlot, args=("SF", conf[0], conf[1], conf[2])) for conf in [("metZAnglePlot", "metDiffPlot",True),("metZAngleUncorrPlot", "metDiffPlot",True),("metZAngleUncorrPlot", "metFracPlot",True), ("metmetAnglePlot", "metDiffPlot",True)]]
    processes += [mp.Process(target=make2DPlot, args=("SF", conf[0], conf[1], conf[2])) for conf in [("metZAngleUncorrPlot_pos_2j", "metDiffPlot",False),("metZAngleUncorrPlot_neg_2j", "metDiffPlot",False),("metZAngleUncorrPlot_pos_3j", "metDiffPlot",False),("metZAngleUncorrPlot_neg_3j", "metDiffPlot",False)]]
    #processes += [mp.Process(target=make2DPlot, args=("SF", conf[0], conf[1], conf[2])) for conf in [("metZAnglePlot", "metDiffPlot",False),("metZAngleUncorrPlot", "metDiffPlot",False),("metZAngleUncorrPlot", "metFracPlot",False), ("metmetAnglePlot", "metDiffPlot",False)]]
    #processes += [mp.Process(target=make2DPlot, args=("SF", conf[0], conf[1], conf[2])) for conf in [("metZAnglePlot", "jzbPlot",False),("metZAngleUncorrPlot", "jzbUncorrPlot",False),("jzbUncorrPlot", "metFracPlot",False), ("jzbUncorrPlot", "metDiffPlot",False),("jzbPlot", "metFracPlot",False), ("jzbPlot", "metDiffPlot",False)]]
    #processes += [mp.Process(target=make2DPlot, args=("SF", conf[0], conf[1], conf[2])) for conf in [("metPlot_3j", "jzbPlot",False),("rawMetPlot_3j", "jzbUncorrPlot",False)]]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    print str(time.time()-t1)+" seconds elapsed"

if __name__ == "__main__":
    main()

