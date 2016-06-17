
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

ERROR_MODE = 1 # 0 only stat, 1 stat and sys, 2 only sys, 3 only pure sys

def produceHistogram(dilepton,mainConfig,stackIt = False):
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
    return fullHist

def getObservedFromSignalRegion(output,region, nJets):
    bkg = getBackgrounds("DY")
    mainConfig = dataMCConfig.dataMCConfig(plot="mllFinePlot",region="Signal%s%dj"%(region,nJets),runName="Run2015_25ns",plotData=False,normalizeToData=False,plotRatio=True,signals=False,useTriggerEmulation=True,personalWork=False,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False)
    
    result = {nJets:{region : {}}}
    
    for dilepton in ["SF", "EE", "MuMu"]:
        HIST_dil = produceHistogram(dilepton,mainConfig,stackIt=False)
        result[nJets][region][dilepton] = HIST_dil
    output.put(result)
    
def getFromSignalRegionNew(output,plotData, region, nJets, dileptons):    
    bkg = getBackgrounds("TT", "DY")

    mainConfig = dataMCConfig.dataMCConfig(plot="mllPlot_neg",region="Signal%s%dj"%(region,nJets),runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=True,signals=False,useTriggerEmulation=True,personalWork=False,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=True, puCorr=True, peakCorr=True)
    mainConfig2 = dataMCConfig.dataMCConfig(plot="mllPlot_neg",region="Signal%s%dj"%(region,nJets),runName="Run2015_25ns",plotData=plotData,normalizeToData=False,plotRatio=True,signals=False,useTriggerEmulation=True,personalWork=False,preliminary=False,forPAS=False,forTWIKI=False,backgrounds=bkg,dontScaleTrig=False,plotSyst=False,doPUWeights=False, responseCorr=True, puCorr=True, peakCorr=True,correctMET=False)
    NEG_OF = produceHistogram("EMu",mainConfig,stackIt=False)
    NEG_OF2 = produceHistogram("EMu",mainConfig2,stackIt=False)
    of_err = ROOT.Double()
    of2_err = ROOT.Double()
    of = NEG_OF.IntegralAndError(1, NEG_OF.GetNbinsX(), of_err)
    of2 = NEG_OF2.IntegralAndError(1, NEG_OF2.GetNbinsX(), of2_err)
    result = {nJets:{region : {"OF":(of+of2, 2*max(float(of2_err),float(of_err)),0)}}}
    
    import corrections
    R = {}  
    R["EE"] = corrections.rEEOF
    R["MuMu"] = corrections.rMMOF
    R["SF"] = corrections.rSFOF
    add = "" if plotData else "MC"
    
    for dilepton in dileptons:
        NEG_dil = produceHistogram(dilepton,mainConfig,stackIt=False)
        NEG2_dil = produceHistogram(dilepton,mainConfig2,stackIt=False)
        dil_err = ROOT.Double()
        dil = NEG_dil.IntegralAndError(1, NEG_dil.GetNbinsX(), dil_err)
        dil2_err = ROOT.Double()
        dil2 = NEG2_dil.IntegralAndError(1, NEG2_dil.GetNbinsX(), dil2_err)
        #print dilepton, region, nJets, dil
        import errorConfig
        sys_err = getattr(errorConfig.DYPredError, region).ERR[nJets]
        fac = float(getattr(getattr(R[dilepton],region.lower()),"val"+add))
        result[nJets][region][dilepton] = (dil+dil2, 2*max(float(dil2_err),float(dil_err)), (dil-fac*of)*sys_err)
    
    output.put(result)

def getErr(tupl):
    tot = 0
    if ERROR_MODE == 0:
        for i in range(1,2):
            tot += (tupl[i])**2
    elif ERROR_MODE == 1:
        for i in range(1,len(tupl)):
            tot += (tupl[i])**2
    elif ERROR_MODE == 2: 
        for i in range(2,len(tupl)):
            tot += (tupl[i])**2
    elif ERROR_MODE == 3: 
        for i in range(len(tupl)-1,len(tupl)):
            tot += (tupl[i])**2
    return tot**0.5

def combineObsMC(dileptons=["SF","EE","MuMu"]):
    import multiprocessing as mp
    output = mp.Queue()
    results = {}
    for nJets in [2,3]:
        results[nJets] = {}
        for region in ["Inclusive", "Forward", "Central"]:
            results[nJets][region] = {}
            for dilepton in dileptons:
                results[nJets][region][dilepton] = None
    
    processes = []
    processes += [mp.Process(target=getObservedFromSignalRegion, args=(output, region, nJets)) for nJets in [2,3] for region in ["Inclusive", "Forward", "Central"]]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    for p in processes:
        out = output.get()
        for key1 in out.keys():
            for key2 in out[key1].keys():
                for key3 in out[key1][key2].keys():
                    results[key1][key2][key3] = out[key1][key2][key3]
                    
    combined = {}
    for region in ["Inclusive", "Forward", "Central"]:
        combined[region] = {}
        for dilepton in dileptons:
            val = results[2][region][dilepton].Clone()
            val.Add(results[3][region][dilepton].Clone(),1)
            combined[region][dilepton] = val
    
    obs_lowMass = {}
    obs_belowZ = {}
    obs_onZ = {}
    obs_aboveZ = {}
    obs_highMass = {}
    
    for dilepton in dileptons:
        obs_lowMass [dilepton] = {}
        obs_belowZ  [dilepton] = {}
        obs_onZ     [dilepton] = {}
        obs_aboveZ  [dilepton] = {}
        obs_highMass[dilepton] = {}
        for region in ["Inclusive", "Forward", "Central"]:
            hist = combined[region][dilepton]
            err_lowMass = ROOT.Double()
            val_lowMass = hist.IntegralAndError(hist.FindBin(20),hist.FindBin(70)-1,err_lowMass)
            err_belowZ = ROOT.Double()
            val_belowZ = hist.IntegralAndError(hist.FindBin(70),hist.FindBin(81)-1,err_belowZ)
            err_onZ = ROOT.Double()
            val_onZ = hist.IntegralAndError(hist.FindBin(81),hist.FindBin(101)-1,err_onZ)
            err_aboveZ = ROOT.Double()
            val_aboveZ = hist.IntegralAndError(hist.FindBin(101),hist.FindBin(120)-1,err_aboveZ)
            err_highMass = ROOT.Double()
            val_highMass = hist.IntegralAndError(hist.FindBin(120),hist.FindBin(300)-1,err_highMass)
            
            obs_lowMass[dilepton][region] = (val_lowMass, float(err_lowMass),0,0)
            obs_belowZ[dilepton][region] = (val_belowZ, float(err_belowZ),0,0)
            obs_onZ[dilepton][region] = (val_onZ, float(err_onZ),0,0)
            obs_aboveZ[dilepton][region] = (val_aboveZ, float(err_aboveZ),0,0)
            obs_highMass[dilepton][region] = (val_highMass, float(err_highMass),0,0)
    
    print "\\begin{tabular}{|l|l|l|l|l|l|l}\\hline"      
    print "%-*s & %-*s& %-*s   & %-*s   & %-*s   & %-*s   & %-*s \\\\\\hline"%(5,"LL", 10, "Region",18,"lowMass",18,"belowZ",18,"zPeak",18,"aboveZ",18,"highMass")
    for dilepton in dileptons:
        for region in ["Inclusive", "Forward", "Central"]:
            print "%-*s &%-*s &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$\\\\\\hline"%(5,dilepton,10,region,8,obs_lowMass[dilepton][region][0],8,getErr(obs_lowMass[dilepton][region]),8,obs_belowZ[dilepton][region][0],8,getErr(obs_belowZ[dilepton][region]),8,obs_onZ[dilepton][region][0],8,getErr(obs_onZ[dilepton][region]),8,obs_aboveZ[dilepton][region][0],8,getErr(obs_aboveZ[dilepton][region]),8,obs_highMass[dilepton][region][0],8,getErr(obs_highMass[dilepton][region]))
    print "\\end{tabular}"
    
    obs = {}
    obs["lowMass"] = obs_lowMass
    obs["belowZ"] = obs_belowZ
    obs["onZ"] = obs_onZ
    obs["aboveZ"] = obs_aboveZ
    obs["highMass"] = obs_highMass
        
    return obs
    
def combineResultsNew(plotData, dileptons=["SF", "EE", "MuMu"]):
    import multiprocessing as mp
    output = mp.Queue()

    results = {}
    for nJets in [2,3]:
        results[nJets] = {}
        for region in ["Inclusive", "Forward", "Central"]:
            results[nJets][region] = {}
            for dilepton in dileptons+["OF"]:
                results[nJets][region][dilepton] = (0,0)

    processes = []
    processes += [mp.Process(target=getFromSignalRegionNew, args=(output,plotData, region, nJets, dileptons)) for nJets in [2,3] for region in ["Inclusive", "Forward", "Central"]]
    
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    
    
    for p in processes:
        out = output.get()
        for key1 in out.keys():
            for key2 in out[key1].keys():
                for key3 in out[key1][key2].keys():
                    results[key1][key2][key3] = out[key1][key2][key3]
    
    combined = {}
    for region in ["Inclusive", "Forward", "Central"]:
        combined[region] = {}
        for dilepton in dileptons+["OF"]:
            val = results[2][region][dilepton][0] + results[3][region][dilepton][0]
            err = ((results[2][region][dilepton][1])**2 + (results[3][region][dilepton][1])**2)**0.5
            sys = (results[2][region][dilepton][2]**2 + results[3][region][dilepton][2]**2)**0.5
            combined[region][dilepton] = (val,err,sys)
    #print combined
    
    import corrections
    R = {}
    outIn = {}
    R["EE"] = corrections.rEEOF
    R["MuMu"] = corrections.rMMOF
    R["SF"] = corrections.rSFOF
    outIn["EE"] = corrections.rOutInEE
    outIn["MuMu"] = corrections.rOutInMM
    outIn["SF"] = corrections.rOutIn
    
    add = "" if plotData else "MC"
    
    lowpred = {}
    belowpred = {}
    zpred = {}
    abovepred = {}
    highpred = {}
    for dilepton in dileptons:
        lowpred[dilepton] = {}
        belowpred[dilepton] = {}
        zpred[dilepton] = {}
        abovepred[dilepton] = {}
        highpred[dilepton] = {}
        for region in ["Inclusive", "Forward", "Central"]:
            fac = getattr(getattr(R[dilepton],region.lower()),"val"+add)
            fac_err = getattr(getattr(R[dilepton],region.lower()),"err"+add)
            val = (combined[region][dilepton][0] - float(fac)*combined[region]["OF"][0])
            val_stat = (combined[region][dilepton][1]**2 + (float(fac)*combined[region]["OF"][1])**2)**0.5
            val_sys = combined[region][dilepton][2]
            val_rsfof = fac_err*combined[region]["OF"][0]
            scale_lowMass = getattr(getattr(outIn[dilepton].lowMass, region.lower()),"val"+add)
            err_lowMass = getattr(getattr(outIn[dilepton].lowMass, region.lower()),"err"+add)
            scale_belowZ = getattr(getattr(outIn[dilepton].belowZ, region.lower()),"val"+add)
            err_belowZ = getattr(getattr(outIn[dilepton].belowZ, region.lower()),"err"+add)
            scale_aboveZ = getattr(getattr(outIn[dilepton].aboveZ, region.lower()),"val"+add)
            err_aboveZ = getattr(getattr(outIn[dilepton].aboveZ, region.lower()),"err"+add)
            scale_highMass = getattr(getattr(outIn[dilepton].highMass, region.lower()),"val"+add)
            err_highMass = getattr(getattr(outIn[dilepton].highMass, region.lower()),"err"+add)
            lowpred[dilepton][region] = (val*scale_lowMass,val_stat*scale_lowMass,val_rsfof*scale_lowMass, val*err_lowMass,val_sys*scale_lowMass)
            belowpred[dilepton][region] = (val*scale_belowZ,val_stat*scale_belowZ,val_rsfof*scale_belowZ, val*err_belowZ,val_sys*scale_belowZ)
            zpred[dilepton][region] = (val,val_stat,val_rsfof,0,val_sys)
            abovepred[dilepton][region] = (val*scale_aboveZ,val_stat*scale_aboveZ,val_rsfof*scale_aboveZ, val*err_aboveZ,val_sys*scale_aboveZ)
            highpred[dilepton][region] = (val*scale_highMass,val_stat*scale_highMass,val_rsfof*scale_highMass, val*err_highMass,val_sys*scale_highMass)
            #print zpred[dilepton][region],lowpred[dilepton][region],highpred[dilepton][region]
    print "\\begin{tabular}{|l|l|l|l|l|l|l}\\hline"
    print "%-*s & %-*s& %-*s   & %-*s   & %-*s   & %-*s   & %-*s \\\\\\hline"%(5,"LL", 10, "Region",18,"lowMass",18,"belowZ",18,"zPeak",18,"aboveZ",18,"highMass")
    for dilepton in dileptons:
        for region in ["Inclusive", "Forward", "Central"]:
            print "%-*s &%-*s &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$ &$%-*.2f\\pm%-*.2f$ \\\\\\hline"%(5,dilepton,10,region,8,lowpred[dilepton][region][0],8,getErr(lowpred[dilepton][region]),8,belowpred[dilepton][region][0],8,getErr(belowpred[dilepton][region]),8,zpred[dilepton][region][0],8,getErr(zpred[dilepton][region]),8,abovepred[dilepton][region][0],8,getErr(abovepred[dilepton][region]),8,highpred[dilepton][region][0],8,getErr(highpred[dilepton][region]))
    print "\\end{tabular}"
    
    pred = {}
    pred["lowMass"] = lowpred
    pred["belowZ"] = belowpred
    pred["onZ"] = zpred
    pred["aboveZ"] = abovepred
    pred["highMass"] = highpred
    
    if dileptons==["SF","EE","MuMu"]:
        from dataMCConfig import dataMCConfig
        import pickle
        name = "shelves/%s_prediction.pkl"%(dataMCConfig.jzbType) if plotData else "shelves/%s_prediction_MC.pkl"%(dataMCConfig.jzbType)
        with open(name, "w") as predictionFile:
            pickle.dump(pred, predictionFile)
            predictionFile.close()
    return pred
    
    
def main():
    import time
    t1 = time.time()
    #dataMCConfig.dataMCConfig.jzbType = "unCorrMet"
    print "PRED DATA"
    pred = combineResultsNew(True)
    #pred = combineResults(True)
    
    print "PRED MC"
    #pred = combineResults(False)
    #pred = combineResultsNew(False)

    print "OBS MC"
    #obs = combineObsMC()


    print str(time.time()-t1)+" seconds elapsed"

if __name__ == "__main__":
    main()
