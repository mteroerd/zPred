
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
        fitFunc.SetRange(max(lowRange,mean-step*sigma),min(mean+step*sigma, upRange))
    return fitFunc    

def produceHistogram(dilepton,mainConfig,stackIt = False,output=None, sort=None):
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
    if mainConfig.plotData:
        fullHist = getDataHist(plot, tree1, tree2)
        fullHist.SetTitle(";%s;%s"%(plot.xaxis,plot.yaxis))
    else:
        counts = {}
        stack = TheStack(processes,mainConfig.runRange.lumi,mainConfig.plot,tree1,tree2,1.0,scaleTree1,scaleTree2,saveIntegrals=True,counts=counts,doTopReweighting=mainConfig.doTopReweighting,theoUncert=mainConfig.theoUncert,doPUWeights=mainConfig.doPUWeights)
        if not stackIt:
            fullHist = stack.theHistogram
            fullHist.SetTitle(";%s;%s"%(plot.xaxis,plot.yaxis))
        else:
            fullHist = stack
            fullHist.theStack.SetTitle(";%s;%s"%(plot.xaxis,plot.yaxis))
            fullHist.theHistogram.SetTitle(";%s;%s"%(plot.xaxis,plot.yaxis))
    if output != None:
        output.put((sort, fullHist))
            
    return fullHist

def doPUCorrection(plotData=True,nJets=2,direction="Central",extraArg=""):
    dilepton = "SF"
    bkg = getBackgrounds("TT","DYAll")

    nPoints = 0
    setTDRStyle()
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetErrorX(0.5)
    
    bins = array('f',[0,7,10,15,31])
    nBins = len(bins)-1
    pointsHist = ROOT.TH1F("h","", nBins,bins)
    
    nVCuts = ["nVert06", "nVert79", "nVert1014", "nVert1530"]
    
    for nVCut in nVCuts:       
        mainConfig = dataMCConfig.dataMCConfig(plot="jzbPlot_puCorr_%dj#jzb50#%s"%(nJets,nVCut),region=direction,runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=False)
        
        histSF = produceHistogram("SF", mainConfig)
        histEM = produceHistogram("EMu", mainConfig)
        if mainConfig.plotData:
            R = getattr(mainConfig.rSFOF, direction.lower()).val
        else:
            R = getattr(mainConfig.rSFOF, direction.lower()).valMC
        histSF.Add(histEM,-R)

        tempHist = histSF.Clone()
        tempHist.Scale(1./tempHist.Integral())
        tempHist.SetMinimum(0)
        can = ROOT.TCanvas()
        tempFunc = ROOT.TF1("temf","gaus",-50,50)
        tempFunc.SetParameter(1,3)
        tempHist.Fit(tempFunc,"RQ")
        tempHist.Draw()
        tempFunc.Draw("same")
        ensurePathExists("fig/%s/PUCorr/%d/%dj/%s/"%(mainConfig.jzbType, mainConfig.correctionMode,nJets,direction))
        can.SaveAs("fig/%s/PUCorr/%d/%dj/%s/%s.pdf"%(mainConfig.jzbType, mainConfig.correctionMode,nJets,direction,nVCut.split("nVert")[1]))
        

        ymean = tempFunc.GetParameter(1)
        yerr = tempFunc.GetParError(1)
        
        nPoints+=1
        pointsHist.SetBinContent(nPoints, ymean)
        pointsHist.SetBinError(nPoints, yerr)
        
     
    template = plotTemplate(mainConfig)
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetErrorX(0.5)
    template.dilepton = dilepton
    template.changeScale = False
    fitLine = TF1("fitLine", "pol1",0,30)
    fitLine.SetParName(1, "#beta")
    fitLine.SetParameter(0,5)
    fitLine.SetParameter(1,0.1)
    pointsHist.Fit(fitLine, "R%s"%(extraArg))
    
    with open("shelves/%s_%d.pkl"%(mainConfig.jzbType,mainConfig.correctionMode), "r+") as corrFile:        
        corrs = pickle.load(corrFile)
        corrFile.seek(0) 
        corrFile.truncate()
        
        corrs[mainConfig.plotData][direction]["pu"][nJets==2] = fitLine.GetParameter(1)
        pickle.dump(corrs, corrFile)
        corrFile.close()
                      
    jetInd = "#geq3j" if nJets == 3 else "=2j"
    template.cutsText = jetInd
    template.labelX = "N_{Vertices}"
    template.labelY = "JZB Peak Position [GeV]"
    template.regionName = direction
                      
    template.setPrimaryPlot(pointsHist, "pe")
    template.addSecondaryPlot(fitLine,"")
                      
    template.draw()   
                      
    template.canvas.Update()
    pointsHist.GetListOfFunctions()
    panel = pointsHist.GetListOfFunctions().FindObject("stats")
    panel.SetX1NDC(0.65)
    panel.SetX2NDC(0.95)
    panel.SetY1NDC(0.14)
    panel.SetY2NDC(0.38)
    template.draw()   
                      
    template.setFolderName("PUCorr/%d"%(mainConfig.correctionMode))
    dataInd = "Data" if plotData else "MC"
    template.saveAs("puSlope_%dj_%s_%s"%(nJets, direction, dataInd))
                      
                      
def main():           
    doPUCorrection(False, 3, "Central")
    doPUCorrection(False, 3, "Forward")
    #for direction in ["Central", "Forward"]:
        #doPUCorrection(True,  2, direction)
        #doPUCorrection(True,  3, direction)
        #doPUCorrection(False, 2, direction)
        #doPUCorrection(False, 3, direction)

if __name__ == "__main__":
    main()
                    
