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
from plotTemplate import plotTemplate
from sampleConfig import getBackgrounds

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

    

def closureTestMC(var, nJets=2, dilepton="SF",region="Inclusive"):
    bkg = getBackgrounds("DY", "DYTauTau")
    mainConfig_a = dataMCConfig.dataMCConfig(plot="%s_pos_%dj"%(var, nJets),region=region,runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=True,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=True, puCorr=True, peakCorr=True,correctMET=False)
    mainConfig_b = dataMCConfig.dataMCConfig(plot="%s_neg_%dj"%(var, nJets),region=region,runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=True,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=True, puCorr=True, peakCorr=True)
    
    template = plotTemplate(mainConfig_a)
    
    POS_SF = produceHistogram(dilepton,mainConfig_a,stackIt=False)
    POS_OF = produceHistogram("EMu",mainConfig_a,stackIt=False)
    NEG_SF = produceHistogram(dilepton,mainConfig_b,stackIt=True)
    NEG_OF = produceHistogram("EMu",mainConfig_b,stackIt=True)
    
    import corrections
    R = {   
            "SF" : getattr(corrections.rSFOF, region.lower()), 
            "EE" : getattr(corrections.rEEOF, region.lower()), 
            "MuMu" : getattr(corrections.rMMOF, region.lower())
        }
    
    
    OBS_SF = POS_SF.Clone()
    OBS_SF.Add(POS_OF,-R[dilepton].valMC)
    PRED_SF = NEG_SF.theHistogram
    PRED_SF.Add(NEG_OF.theHistogram,-R[dilepton].valMC)

    pred_err = ROOT.Double()
    pred = PRED_SF.IntegralAndError(0,PRED_SF.GetNbinsX(), pred_err)
    obs_err = ROOT.Double()
    obs = OBS_SF.IntegralAndError(0,OBS_SF.GetNbinsX(), obs_err)
    
    PRED_SF.SetLineColor(ROOT.kBlack)
    PRED_SF.SetLineWidth(2)
    PRED_SF.SetMarkerSize(0)

    template.setPrimaryPlot(OBS_SF, "PE")
    template.maximumScale = 100
    template.logY = True
    template.addSecondaryPlot(PRED_SF, "HISTE")
    
    template.regionName = region+", Z+jets only"
    jetInd = "#geq3j" if (nJets == 3) else "=2j" 
    if "ATLAS" in var:
        jetInd = "#geq3j, ATLAS" if (nJets == 3) else "=2j, ATLAS" 
    template.cutsText = jetInd
    template.cutsSize/=0.7

    template.dilepton = dilepton
    template.regionSize/=0.7
    template.ratioLabel = "JZB>0/JZB<0"
    template.hasRatio = True
    
    if "met" in var:
        import errorConfig
        ERR = getattr(errorConfig.DYPredError, region).ERR[nJets]
        template.addRatioErrorBySize("Mismatch of spectra", ERR, ROOT.kGreen-6, 1001, False)
    elif "jzb" in var and nJets==3:
        ERR = 0.3
        template.addRatioErrorBySize("Mismatch of spectra", ERR, ROOT.kGreen-6, 1001, False)
    elif "jzb" in var and nJets==2:
        ERR = 0.5
        template.addRatioErrorBySize("Mismatch of spectra", ERR, ROOT.kGreen-6, 1001, False)
    
    template.draw()
    template.plotPad.cd()
            
    leg = TLegend(0.55,0.75,0.9,0.9)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.045)
    leg.AddEntry(OBS_SF,   "JZB>0", "PE")
    leg.AddEntry(PRED_SF,  "JZB<0","l")
    leg.Draw("same")
    
    if "met" in var:
        cut = 100 if nJets==3 else 150
        cutL = ROOT.TLine(cut,0,cut,1e2)
        cutL.SetLineStyle(2)
        cutL.SetLineWidth(2)
        cutL.Draw("same")
        
        
    template.setFolderName("Tests/%d"%(mainConfig_a.correctionMode))
    template.saveAs("%s_%dj_%s_%s_Closure_Test_MC"%(var,nJets,dilepton,region))
   
def compareJZBforSamples():
    bkg = getBackgrounds("DY")
    mainConfigDY = dataMCConfig.dataMCConfig(plot="jzbPlot_3j_allMasses",region="Inclusive",runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=True, puCorr=True, peakCorr=True,correctMET=False)
    
    bkg = getBackgrounds("TT")
    mainConfigTT = dataMCConfig.dataMCConfig(plot="jzbPlot_3j_allMasses",region="Inclusive",runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=False,signals=False,useTriggerEmulation=True,personalWork=True,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=True, puCorr=True, peakCorr=True,correctMET=False)

    SUSYTreesEE = readTrees("/user/teroerde/trees/sw7414v2005/", "EE")
    SUSYTreesMuMu = readTrees("/user/teroerde/trees/sw7414v2005/", "MuMu")
    
    EE = SUSYTreesEE["T6bbllslepton_msbottom_400_mneutralino_150"]
    MuMu = SUSYTreesMuMu["T6bbllslepton_msbottom_400_mneutralino_150"]
    
    template = plotTemplate(mainConfigTT)
    
    histDY = produceHistogram("SF",mainConfigDY)
    histTTBAR = produceHistogram("SF",mainConfigTT)
    
    histEE = createHistoFromTree(EE, mainConfigTT.plot.variable, mainConfigTT.plot.cuts,80,-200,200, binning=[])
    histMuMu = createHistoFromTree(MuMu, mainConfigTT.plot.variable, mainConfigTT.plot.cuts,80,-200,200, binning=[])
    histSUSY = histEE.Clone()
    histSUSY.Add(histMuMu.Clone(),1)
    
    #return
    histDY.SetLineColor(401)
    histTTBAR.SetLineColor(855)
    histSUSY.SetLineColor(ROOT.kRed)
    
    histDY.SetLineWidth(2)
    histTTBAR.SetLineWidth(2)
    histSUSY.SetLineWidth(2)
    
    histDY.Scale(1./histDY.Integral())
    histTTBAR.Scale(1./histTTBAR.Integral())
    histSUSY.Scale(1./histSUSY.Integral())
    
    histDY.SetMinimum(1e-3)
    
    template.setPrimaryPlot(histDY,"hist")
    template.addSecondaryPlot(histTTBAR, "hist")
    template.addSecondaryPlot(histSUSY, "hist")
    template.lumiInt = None

    template.dilepton = "SF"
    template.logY = True
    template.maximumScale = 0.5e2
    template.regionName = "Inclusive"
    template.cutsText = "#geq3j" 
    
    template.draw()

    leg = ROOT.TLegend(0.55,0.73,0.93,0.84)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(histDY, "LO DY+jets", "l")
    leg.AddEntry(histTTBAR, "t#bar{t} Powheg", "l")
    leg.AddEntry(histSUSY, "SUSY Signal", "l")
    leg.Draw("same")
    

    template.setFolderName("JZBCompare")
    template.saveAs("comparison")
    
def main():
    import multiprocessing as mp
    import time
    t1 = time.time()
    processes = []
    #processes += [mp.Process(target=closureTestMC, args=(var, nJets)) for nJets in [2,3] for var in ["metDiffPlot"]]
    #processes += [mp.Process(target=closureTestMC, args=(var, nJets,dil,region)) for dil in ["SF"] for nJets in [2,3] for var in ["jzbPlot"] for region in ["Forward", "Central"]]
    #processes += [mp.Process(target=closureTestMC, args=(var, nJets,dil,region)) for dil in ["SF"] for nJets in [2,3] for var in ["metPlot"] for region in ["Forward", "Central"]]
    processes += [mp.Process(target=closureTestMC, args=(var, nJets,dil,region)) for dil in ["SF"] for nJets in [2,3] for var in ["jzbPlot_ATLAS"] for region in ["Forward", "Central"]]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    #compareJZBforSamples()
    print str(time.time()-t1)+" seconds elapsed"

if __name__ == "__main__":
    main()
