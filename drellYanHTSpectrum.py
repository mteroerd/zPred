
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


def produceHistogram(dilepton,mainConfig,stackIt = False,output=None, sort=None):
    eventCounts = totalNumberOfGeneratedEvents(mainConfig.dataSetPath)  

    processes = []
    for background in mainConfig.backgrounds:
        processes.append(Process(getattr(Backgrounds,background),eventCounts))

        #~ if mainConfig.plotSignal:
            #~ processes.append(Signal)
    nEvents=-1
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
    
    #plotPad.SetLogy()
    if mainConfig.useTriggerEmulation or mainConfig.DontScaleTrig:
        scaleTree2 = 1.0
        scaleTree1 = 1.0                
    if mainConfig.plotData:
        fullHist = getDataHist(plot, tree1, tree2)
        fullHist.SetTitle(";%s;%s"%(plot.xaxis,plot.yaxis))
        #fullHist.Draw("PE")
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
    #hCanvas.Print("fig/Closure/"+mainConfig.plot.filename%("_"+nameModifier),)
 

def plotHTSpectrum(nJets, dilepton="SF"):
    bkg = ["DrellYanLO", "DrellYanLOHT0to100"]
    if nJets==0:
        mainConfig = dataMCConfig.dataMCConfig(plot="genHTPlotFine",region="Inclusive",runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False)
    else:
        mainConfig = dataMCConfig.dataMCConfig(plot="genHTPlotFine_%dj"%(nJets),region="Inclusive",runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False)
    
    from plotTemplate import plotTemplate
    template = plotTemplate(mainConfig)
    template.logY = True
    template.maximumScale = 1e1
    
    SF = produceHistogram(dilepton,mainConfig,stackIt=False)
    SF.SetLineColor(ROOT.kBlack)
    SF.SetMarkerSize(0)

    template.setPrimaryPlot(SF, "HISTE")
    
    if nJets==0:
        jetind= ""
    elif nJets==3:
        jetind = "#geq3j"
    elif nJets==2:
         jetind = "=2j"
    template.cutsText = "#splitline{%s}{%s}"%("Z+jets only",jetind)
    template.cutsPosX = 0.80
    template.cutsPosY = 0.87
    
    
    template.draw()
    
    arrs = []
    for pos in [100,200,400,600]:
        arr = ROOT.TArrow(pos,SF.GetBinContent(SF.FindBin(pos)-1)*4,pos,SF.GetBinContent(SF.FindBin(pos)-1)*1.4)
        arrs.append(arr)
        arr.SetNDC(False)
        arr.SetArrowSize(0.03)
        arr.SetAngle(40)
        arr.Draw()
    
    template.setFolderName("DrellYanHTSpectrum")
    template.saveAs("spectrum_%dj_%s"%(nJets,dilepton))

def main():
    import multiprocessing as mp
    import time
    t1 = time.time()
    
    plotHTSpectrum(0,"SF")
    plotHTSpectrum(2,"SF")
    plotHTSpectrum(3,"SF")
   
    print str(time.time()-t1)+" seconds elapsed"

if __name__ == "__main__":
    main()
