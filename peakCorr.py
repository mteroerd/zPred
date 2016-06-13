
import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np


import pickle

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
from plotTemplate import plotTemplate, plotTemplate2D
from sampleConfig import getBackgrounds

def gaussianFit(histo, lowRange, upRange,step):
    fitFunc = TF1("name", "gaus", lowRange, upRange)
    fitFunc.SetParName(1, "#mu")
    #first guess
    mean = histo.GetBinCenter(histo.GetMaximumBin())
    #mean = histo.GetMean()
    #second guess
    sigma = (upRange-lowRange)/4
    histo.Fit(fitFunc,"QLLN","", mean-sigma, mean+sigma)
    mean = fitFunc.GetParameter(1)
    sigma = fitFunc.GetParameter(2)
    for i in range(5):
        histo.Fit(fitFunc, "QLL0","", mean-step*sigma, mean+step*sigma)
        mean = fitFunc.GetParameter(1)
        sigma = fitFunc.GetParameter(2)
        fitFunc.SetRange(max(lowRange+0.001,mean-step*sigma),min(mean+step*sigma, upRange-0.001))
    return fitFunc    

def doPeakCorrection(plotData=True, nJets=2, direction="Central", extraArg=""):
    dilepton = "SF"
    bkg = getBackgrounds("TT", "DY", "DYTauTau")

    mainConfig = dataMCConfig.dataMCConfig(plot="jzbPlot_peakCorr_%dj#jzb50"%(nJets),region=direction,runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False,responseCorr=True, puCorr=True)

    eventCounts = totalNumberOfGeneratedEvents(mainConfig.dataSetPath)  

    processes = []
    for background in mainConfig.backgrounds:
        processes.append(Process(getattr(Backgrounds,background),eventCounts))

    treeEE = readTrees(mainConfig.dataSetPath, "EE")
    treeMuMu = readTrees(mainConfig.dataSetPath, "MuMu")
    treeEMu = readTrees(mainConfig.dataSetPath, "EMu")
 
    mainConfig.plot.addDilepton(dilepton)    
    plot = mainConfig.plot
    
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
    
    template = plotTemplate(mainConfig)
    template.maximumScale = 1.1
    template.dilepton = dilepton
    template.regionName = direction
    gStyle.SetOptFit(111)  
    
    if mainConfig.plotData:
        fullHist = getDataHist(plot, tree1, tree2)
        template.setPrimaryPlot(fullHist, "PE")
    else:
        counts = {}
        stack = TheStack(processes,mainConfig.runRange.lumi,mainConfig.plot,tree1,tree2,1.0,scaleTree1,scaleTree2,saveIntegrals=True,counts=counts,doTopReweighting=mainConfig.doTopReweighting,theoUncert=mainConfig.theoUncert,doPUWeights=mainConfig.doPUWeights)
        fullHist = stack.theHistogram
        fullHist.SetFillStyle(0)
        fullHist.SetMarkerSize(0)
        template.setPrimaryPlot(fullHist, "HISTE")
    
    if mainConfig.plotData:
        eMuHist = getDataHist(plot, treeEMu, "None")
    else:
        counts = {}
        eMustack = TheStack(processes,mainConfig.runRange.lumi,mainConfig.plot,treeEMu,"None",1.0,scaleTree1,scaleTree2,saveIntegrals=True,counts=counts,doTopReweighting=mainConfig.doTopReweighting,theoUncert=mainConfig.theoUncert,doPUWeights=mainConfig.doPUWeights)
        eMuHist = eMustack.theHistogram
    
    if mainConfig.plotData:
        R = getattr(mainConfig.rSFOF, direction.lower()).val
    else:
        R = getattr(mainConfig.rSFOF, direction.lower()).valMC
    fullHist.Add(eMuHist,-R)
    
    #fitFunc = gaussianFit(fullHist,-50, 50, 1.6)
    fitFunc = ROOT.TF1("fitf", "gaus",-50,50)
    fullHist.Fit(fitFunc,"RQ")
    template.addSecondaryPlot(fitFunc)
    
    fileName = "%s_%d.pkl"%(mainConfig.jzbType,mainConfig.correctionMode)
    with open("shelves/"+fileName, "r+") as corrFile:
        corrs = pickle.load(corrFile)
        corrFile.seek(0)
        corrFile.truncate()
        
        corrs[mainConfig.plotData][direction]["peak"][nJets==2] = fitFunc.GetParameter(1)
        pickle.dump(corrs, corrFile)
        corrFile.close()
    
    jetInd = "#geq3j" if nJets == 3 else "=2j"
    template.cutsText = jetInd
    
    template.draw()
    
    template.canvas.Update()
    fullHist.GetListOfFunctions()
    panel = fullHist.GetListOfFunctions().FindObject("stats")
    panel.SetY1NDC(0.5)
    panel.SetY2NDC(0.65)
    panel.SetX1NDC(0.7)
    panel.SetX2NDC(1-template.marginRight-0.01)
    template.draw()
    indicator = "Data" if mainConfig.plotData else "MC"
    
    template.setFolderName("PeakCorr/%d"%(mainConfig.correctionMode))
    template.saveAs("peakCorr_%dj_%s_%s"%(nJets,direction, indicator))
    

def main():
    for direction in ["Central", "Forward"]:
        doPeakCorrection(True,  2, direction)
        doPeakCorrection(True,  3, direction)
        doPeakCorrection(False, 2, direction)
        doPeakCorrection(False, 3, direction)

if __name__ == "__main__":
    main()
