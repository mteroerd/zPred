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
from plotTemplate import plotTemplate2D
from sampleConfig import getBackgrounds

def produceHistogram2D(dilepton,mainConfig):
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
    
    #plotPad.SetLogy()
    if mainConfig.useTriggerEmulation or mainConfig.DontScaleTrig:
        scaleTree2 = 1.0
        scaleTree1 = 1.0                
    if mainConfig.plotData:
        fullHist = getData2DHist(plot, plot2, tree1, tree2)
    else:
        MCHists = []
        for process in processes:
            MCHists.append(process.createCombined2DHistogram(mainConfig.runRange.lumi,plot,plot2,tree1,tree2,1,scaleTree1,scaleTree2,doTopReweighting=mainConfig.doTopReweighting, doPUWeights=mainConfig.doPUWeights))
        fullHist = MCHists[0].Clone()
        fullHist.Reset()
        for hist in MCHists:
            fullHist.Add(hist,1)
    fullHist.SetTitle(";%s;%s"%(plot.xaxis,plot2.xaxis))
    return fullHist

def METJZBPlot(nJets, responseCorr, puCorr, peakCorr, metCorr):
    bkg = getBackgrounds("DY")

    mainConfig_pos = dataMCConfig.dataMCConfig(plot="met2DPlot_%dj"%(nJets), plot2="jzb2DPlot_pos",region="Inclusive",runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=responseCorr, puCorr=puCorr, peakCorr=peakCorr, correctMET=False)
    mainConfig_neg = dataMCConfig.dataMCConfig(plot="met2DPlot_%dj"%(nJets), plot2="jzb2DPlot_neg",region="Inclusive",runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=responseCorr, puCorr=puCorr, peakCorr=peakCorr, correctMET=metCorr)
    
    template = plotTemplate2D(mainConfig_pos)
    template.labelZ = "Events"
    template.logZ = True
    
    hist_pos = produceHistogram2D("SF", mainConfig_pos)
    hist_neg = produceHistogram2D("SF", mainConfig_neg)
    
    hist_pos.Add(hist_neg,1)

    template.setPrimaryPlot(hist_pos, "colz")

    if True:
        prof = hist_pos.ProfileX()
        template.addSecondaryPlot(prof)
    
    line = ROOT.TLine(0,0,50,0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.SetNDC(False)
    template.addSecondaryPlot(line)

    jetInd = "#geq3j" if (nJets == 3) else "=2j" 
    template.regionName = "Z+Jets only"
    template.cutsText = jetInd
    template.dilepton = "SF"
    
    responseText = "Response corrections: No"
    if responseCorr:
            responseText = "Response corrections: Yes"
    puText = "Pile-Up corrections: No"
    if puCorr:
            puText = "Pile-Up corrections: Yes"
    peakText = "Peak corrections: No"
    if peakCorr:
            peakText = "Peak corrections: Yes"
    metText = "Propagation to E_{T}^{miss}: No"
    if metCorr:
            metText = "Propagation to E_{T}^{miss}: Yes"
            
    
    responseLatex = ROOT.TLatex()
    responseLatex.SetNDC(True)
    responseLatex.SetTextFont(42)
    responseLatex.SetTextSize(0.035)
    responseLatex.SetText(0.17, 0.3, responseText)
    puLatex = responseLatex.Clone()
    puLatex.SetText(0.17, 0.26, puText)
    peakLatex = responseLatex.Clone()
    peakLatex.SetText(0.17, 0.22, peakText)
    metLatex = responseLatex.Clone()
    metLatex.SetText(0.17, 0.18, metText)
    
    template.addSecondaryPlot(responseLatex)
    template.addSecondaryPlot(puLatex)
    template.addSecondaryPlot(peakLatex)
    template.addSecondaryPlot(metLatex)
    
    nameModifier = "MC"
    
    template.draw()
    
    template.setFolderName("metJZBTests")
    template.saveAs("response%s_pileup%s_peak%s_met%s_%dj_%s"%(responseCorr,puCorr,peakCorr, metCorr, nJets, nameModifier))

def main():
    import multiprocessing as mp
    import time
    t1 = time.time()

    processes = []
    for nJets in [2,3]:
        processes.append(mp.Process(target=METJZBPlot, args =(nJets, False, False, False, False )))
        processes.append(mp.Process(target=METJZBPlot, args =(nJets, True,  False, False, False )))
        processes.append(mp.Process(target=METJZBPlot, args =(nJets, True,  True,  False, False )))
        processes.append(mp.Process(target=METJZBPlot, args =(nJets, True,  True,  True,  False )))
        processes.append(mp.Process(target=METJZBPlot, args =(nJets, True,  True,  True,  True  )))

    for p in processes:
        p.start()
    for p in processes:
        p.join()
        
    print str(time.time()-t1)+" seconds elapsed"

if __name__ == "__main__":
    main()

