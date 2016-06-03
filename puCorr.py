
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
#gROOT.SetBatch(True)
from helpers import *   
import math

from centralConfig import plotLists
from plotTemplate import plotTemplate, plotTemplate2D
from sampleConfig import getBackgrounds

def gaussianFit(histo, lowRange, upRange,step):
    fitFunc = TF1("name", "gaus", lowRange, upRange)
    fitFunc.SetParName(1, "#mu")
    #first guess
    mean = histo.GetBinCenter(histo.FindBin(histo.GetMean()))
    fitFunc.SetParameter(1,mean)
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

def doPUCorrection(plotData=True,nJets=2,direction="Central",extraArg=""):
    dilepton = "SF"
    bkg = getBackgrounds("TT", "DY")
    mainConfig = dataMCConfig.dataMCConfig(plot="nVerticesPlot",plot2="jzbPlot_puCorr_%dj"%(nJets),region=direction,runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=False)
    
    eventCounts = totalNumberOfGeneratedEvents(mainConfig.dataSetPath)  
    processes = []
    for background in mainConfig.backgrounds:
        processes.append(Process(getattr(Backgrounds,background),eventCounts))
           
    ROOT.gStyle.SetOptStat(0)
    
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
    
    template = plotTemplate2D(mainConfig)
    template.dilepton = dilepton
    template.logZ = True
    template.labelZ = "Events"
    jetInd = "#geq3j" if nJets == 3 else "=2j"
    template.cutsText = jetInd
    
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
    
    if mainConfig.plotData:
        eMuHist = getData2DHist(plot,plot2, treeEMu, "None")
    else:
        eMuMCHists = []
        for process in processes:
            eMuMCHists.append(process.createCombined2DHistogram(mainConfig.runRange.lumi,plot,plot2,treeEMu,"None",1,scaleTree1,scaleTree2,doTopReweighting=mainConfig.doTopReweighting, doPUWeights=mainConfig.doPUWeights))
        eMuHist = eMuMCHists[0].Clone()
        eMuHist.Reset()
        for hist in eMuMCHists:
            eMuHist.Add(hist,1)
    
    if mainConfig.plotData:
        R = getattr(mainConfig.rSFOF, direction.lower()).val
    else:
        R = getattr(mainConfig.rSFOF, direction.lower()).valMC
    fullHist.Add(eMuHist,-R)
    
    template.setPrimaryPlot(fullHist, "COLZ")
    
    template.draw()
    template.setFolderName("PUCorr/%d"%(mainConfig.correctionMode))
    dataInd = "Data" if plotData else "MC"
    template.saveAs("pu2Dplot_%dj_%s_%s"%(nJets,direction,dataInd))
    
    
    template = plotTemplate(mainConfig)
    template.dilepton = dilepton
    template.changeScale = False
    
    errPoints = TGraphAsymmErrors()
    nBins = fullHist.GetYaxis().GetNbins()
    #binning = [-1,5,8,11,14,17,30]
    binning = [-1,6,9,14,30]
    
    nPoints = 0
    ROOT.gStyle.SetOptFit(111)
    for i in range(0, len(binning)-1):       
        tempHist = fullHist.ProjectionY("", fullHist.GetXaxis().FindBin(binning[i]+1),fullHist.GetXaxis().FindBin(binning[i+1]), "e")
        #if tempHist.Integral() < 50:
            #tempHist.Rebin(4)
        #tempFunc = gaussianFit(tempHist,-40,40, 2)
        tempFunc = ROOT.TF1("temf","gaus",-50,50)
        tempHist.Fit(tempFunc,"R")
        ymean = tempFunc.GetParameter(1)
        yerr = tempFunc.GetParError(1)
        fullHist.GetXaxis().SetRange(fullHist.GetXaxis().FindBin(binning[i]+1),fullHist.GetXaxis().FindBin(binning[i+1]))
        xmean = fullHist.GetMean(1)
        errPoints.SetPoint(nPoints, xmean, ymean)
        errPoints.SetPointError(nPoints, xmean-binning[i]-0.5, binning[i+1]-xmean+0.5,yerr,yerr)
        nPoints+=1
        
    fullHist.GetXaxis().SetRange(0,30)
    fitLine = TF1("fitLine", "pol1",0,30)
    fitLine.SetParName(1, "#beta")
    fitLine.SetParameter(0,0)
    fitLine.SetParameter(1,0.5)
    errPoints.Fit(fitLine, "R%s"%(extraArg))
    
    with open("shelves/%s_%d.pkl"%(mainConfig.jzbType,mainConfig.correctionMode), "r+") as corrFile:        
        corrs = pickle.load(corrFile)
        corrFile.seek(0) 
        corrFile.truncate()
        
        corrs[mainConfig.plotData][direction]["pu"][nJets==2] = fitLine.GetParameter(1)
        pickle.dump(corrs, corrFile)
        corrFile.close()
    
    jetInd = "#geq3j" if nJets == 3 else "=2j"
    template.cutsText = jetInd
    template.labelX = fullHist.GetXaxis().GetTitle()
    template.labelY = "JZB Peak Position [GeV]"
    template.regionName = direction

    template.setPrimaryPlot(errPoints, "APE")
    template.addSecondaryPlot(fitLine,"")
    
    template.draw()
    
    template.canvas.Update()
    errPoints.GetListOfFunctions()
    panel = errPoints.GetListOfFunctions().FindObject("stats")
    panel.SetX1NDC(0.65)
    panel.SetX2NDC(0.95)
    panel.SetY1NDC(0.14)
    panel.SetY2NDC(0.38)
    template.draw()
    
    template.setFolderName("PUCorr/%d"%(mainConfig.correctionMode))
    dataInd = "Data" if plotData else "MC"
    template.saveAs("puSlope_%dj_%s_%s"%(nJets, direction, dataInd))
    

def main():
    #doPUCorrection(True, 3, "Forward")
    for direction in ["Central", "Forward"]:
        doPUCorrection(True,  2, direction)
        doPUCorrection(True,  3, direction)
        doPUCorrection(False, 2, direction)
        doPUCorrection(False, 3, direction)

if __name__ == "__main__":
    main()
