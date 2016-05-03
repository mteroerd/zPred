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

from ROOT import TCanvas, TPad, TH1F, TH1I, THStack, TLegend, TMath, gROOT, TF1, TGraphErrors,gStyle
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
from plotTemplate import plotTemplate

def jzbBinned(dilepton, plotData=True, corrs=False):
    bkg = ["TT_Powheg","DrellYanLO", "DrellYanLOHT0to100"]
    mainConfig = dataMCConfig.dataMCConfig(plot="jzbPlot",region="Inclusive",runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True, preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False)

    eventCounts = totalNumberOfGeneratedEvents(mainConfig.dataSetPath)  
    processes = []
    for background in mainConfig.backgrounds:
        processes.append(Process(getattr(Backgrounds,background),eventCounts))

    
    ROOT.gStyle.SetOptStat(0)

    treeEE = readTrees(mainConfig.dataSetPath, "EE")
    treeMuMu = readTrees(mainConfig.dataSetPath, "MuMu")
    treeEMu = readTrees(mainConfig.dataSetPath, "EMu")
    
    template = plotTemplate(mainConfig)
    template.labelY = "normalized entries"
    
    plotList = ["jzbPlotPt100","jzbPlotPt150","jzbPlotPt200","jzbPlotPtinf"]
    cutLabelList = []
    histList = []
    colorList = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]
    for toPlot in plotList:
        mainConfig = dataMCConfig.dataMCConfig(plot=toPlot,region="Inclusive",runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=False,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=mainConfig.doPUWeights, responseCorr=corrs)
        mainConfig.plot.addDilepton(dilepton)    
        plot = mainConfig.plot
        cutLabelList.append(plot.label3.split('GeV')[0]+"GeV")
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

            histList.append(fullHist)
        else:
            MCHists = []
            for process in processes:
                MCHists.append(process.createCombinedHistogram(mainConfig.runRange.lumi,plot,tree1,tree2,1,scaleTree1,scaleTree2, doTopReweighting=mainConfig.doTopReweighting ,doPUWeights=mainConfig.doPUWeights))
            fullHist = MCHists[0].Clone()
            fullHist.Reset()
            for hist in MCHists:
                fullHist.Add(hist,1)
            histList.append(fullHist)
        
        fullHist.Scale(1./fullHist.Integral(0,fullHist.GetNbinsX()+1))
        
        
    leg = TLegend(0.7,0.7,0.95,0.95)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)    
        
    for i in range(0, len(histList)):
        histList[i].SetLineColor(colorList[i])
        histList[i].SetFillStyle(0)
        histList[i].GetYaxis().SetTitleOffset(1.2)
        leg.AddEntry(histList[i],cutLabelList[i],"l")
        if i == 0:
            template.setPrimaryPlot(histList[i],"hist")
            histList[i].GetXaxis().SetRangeUser(-65,65)
            
        else:
            template.addSecondaryPlot(histList[i],"hist")
    
    template.draw()
    
    leg.Draw("same")
    
    corrString = "corr" if corrs else "uncorr"
    dataInd = "Data" if plotData else "MC"
    
    template.setFolderName("JZBBinned")
    template.saveAs("binned_%s_%s"%(corrString, dataInd))
    
from plotTemplate import plotTemplate2D    
def doResponseCorrection(dilepton,plotData,extraArg=""):
    bkg = ["TT_Powheg","DrellYanLO", "DrellYanLOHT0to100"]
    mainConfig = dataMCConfig.dataMCConfig(plot2="responsePlot",plot="ptllresponsePlot",region="Inclusive",runName="Run2015_25ns",plotData=plotData,personalWork=True,backgrounds=bkg)

    template = plotTemplate2D(mainConfig)        
    template.dilepton = dilepton
    template.logZ = True
    template.redrawPrimary = False
    
    template.labelZ = "Events"
    
    eventCounts = totalNumberOfGeneratedEvents(mainConfig.dataSetPath)  
    processes = []
    for background in mainConfig.backgrounds:
        processes.append(Process(getattr(Backgrounds,background),eventCounts))

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
    
    errPoints = fullHist.ProfileX("", 2, fullHist.GetNbinsY()-1)
    fitFunc = TF1("fitfunc", "pol0",100,400)
    errPoints.Fit(fitFunc, "R%s"%(extraArg))
    
    template.setPrimaryPlot(fullHist, "COLZ")
    template.addSecondaryPlot(errPoints, "PE")
    template.addSecondaryPlot(fitFunc,"")
    
    template.draw()
    
    with open("shelves/%s.pkl"%(mainConfig.jzbType), "r+") as corrFile:
        corrs = pickle.load(corrFile)
        
        corrFile.seek(0)
        corrFile.truncate()
        corrs[mainConfig.plotData]["response"] = 1-fitFunc.GetParameter(0)
        pickle.dump(corrs, corrFile)
        corrFile.close()
    
    leg = TLegend(0.79-template.marginRight,0.65,0.99-template.marginRight,0.75)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)        
    leg.AddEntry(errPoints,"Mean R_{JZB} Value","PE")
    leg.AddEntry(fitFunc, "Constant Fit:","l")
    leg.AddEntry(fitFunc,"R_{JZB}=%.4f"%(fitFunc.GetParameter(0)),"")
    leg.Draw("same")
    
    if mainConfig.plotData:
        name = "Response_Data"
    else:
        name = "Response_MC"
    
    template.setFolderName("Response")
    template.saveAs(name)

def main():
    import multiprocessing as mp
    import time
    t1 = time.time()

    processes = []
    #processes += [mp.Process(target=jzbBinned, args=("SF", plotData, responseCorr)) for plotData in [False, True] for responseCorr in [False, True]]
    processes += [mp.Process(target=doResponseCorrection, args=("SF", plotData)) for plotData in [False, True]]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    print str(time.time()-t1)+" seconds elapsed"
    

if __name__ == "__main__":
    main()