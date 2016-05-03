import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np


from ROOT import TMath
from locations import locations
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

from centralConfig import regionsToUse

def produceHistograms(plot):
    dataSetPath = locations.dataSetPath
    histEE      = getDataHist(plot, readTrees(dataSetPath, "EE")  )
    histMuMu    = getDataHist(plot, readTrees(dataSetPath, "MuMu"))
    histEMu     = getDataHist(plot, readTrees(dataSetPath, "EMu") )
    return histEE, histMuMu, histEMu

def getLines(yMin,yMax, xPos = [70.,81., 101., 120]):
    from ROOT import TLine, kGray
    result = []
    for x in xPos:
        result.append(TLine(x, yMin, x,yMax))
        result[-1].SetLineWidth(1)
        result[-1].SetLineColor(kGray+2)
        result[-1].SetLineStyle(2)
    return result

from plotTemplate import plotTemplate
def makePrediction(region, dil, plot, TT_full, R, data, histTT, histDY, histDYScale):
    template = plotTemplate(TT_full)
    
    ####################################################################
    # ttbar
    ####################################################################
    
    TTpred = histTT.Clone()
    TTpredSys = TTpred.Clone()
    
    #scale OF up with RSFOF
    TTpred.Scale(R[0])
    for i in range(1, TTpred.GetNbinsX()+1):
            TTpredSys.SetBinContent(i,0)
            TTpredSys.SetBinError(i,R[1]*histTT.GetBinContent(i))
     
    ####################################################################
    # z+jets
    ####################################################################
    
    import pickle
    zpred = None
    with open("shelves/%s_prediction.pkl"%(dataMCConfig.dataMCConfig.jzbType),"r") as zpredFile:
        zpred = pickle.load(zpredFile)
        zpredFile.close()
    
    #factors for relative statistical and systematic errors
    zpred_rel_stat = zpred["onZ"][dil][region][1] / zpred["onZ"][dil][region][0]
    zpred_rel_sys = sum(zpred["onZ"][dil][region][2:]) / zpred["onZ"][dil][region][0]
    
    DYpred = histDY.Clone()
    DYpredSys = histDY.Clone()
    
    #scale histogram to match number of predicted events in onZ region
    DYpred.Scale(zpred["onZ"][dil][region][0]/histDYScale.Integral(histDYScale.FindBin(81), histDYScale.FindBin(101)))
    
    from centralConfig import mllBins
    from corrections import rOutIn
    
    
    #onZ errors
    for i in range(DYpredSys.FindBin(mllBins.onZ.low+0.01), DYpredSys.FindBin(mllBins.onZ.high-0.01)):
        DYpred.SetBinError(i,zpred_rel_stat*DYpred.GetBinContent(i))
        DYpredSys.SetBinContent(i,0)
        DYpredSys.SetBinError(i,zpred_rel_sys*DYpred.GetBinContent(i))
        
    #offZ errors
    offZRegions = [(mllBins.lowMass.low,mllBins.lowMass.high, "lowMass"), (mllBins.belowZ.low,mllBins.belowZ.high, "belowZ"),(mllBins.aboveZ.low,mllBins.aboveZ.high, "aboveZ"),(mllBins.highMass.low,plot["DY"].lastBin, "highMass")]
    for bounds in offZRegions:
        for i in range(DYpredSys.FindBin(bounds[0]),DYpredSys.FindBin(bounds[1])):
            DYpred.SetBinError(i,(zpred_rel_stat)*DYpred.GetBinContent(i))
            DYpredSys.SetBinContent(i,0) 
            DYpredSys.SetBinError(i,(zpred_rel_sys**2 + (getattr(getattr(rOutIn, bounds[2]), region.lower()).err/getattr(getattr(rOutIn, bounds[2]), region.lower()).val)**2)**(0.5) *DYpred.GetBinContent(i))
    
    ####################################################################
    #combine DY and TT prediction
    ####################################################################
    
    pred  = ROOT.TH1F("pred%s%s"%(region,dil),"",plot["DY"].nBins,plot["DY"].firstBin, plot["DY"].lastBin)
    predErr = ROOT.TGraphAsymmErrors()
    predUp   = pred.Clone()
    predDown = pred.Clone()
    
    #calculate errors
    for i in range(1, pred.GetNbinsX()+1):
        pred.SetBinContent(i,TTpred.GetBinContent(i)+DYpred.GetBinContent(i))
        pred.SetBinError(i, (TTpred.GetBinError(i)**2 + DYpred.GetBinError(i)**2)**0.5)
        errFull = (TTpred.GetBinError(i)**2 + TTpredSys.GetBinError(i)**2 + DYpred.GetBinError(i)**2 + DYpredSys.GetBinError(i)**2)**0.5
        errSyst = (TTpredSys.GetBinError(i)**2 + DYpredSys.GetBinError(i)**2)**0.5
        
        #total uncertainty graph
        predErr.SetPoint(i, plot["DY"].firstBin + (i-0.5)*((plot["DY"].lastBin-plot["DY"].firstBin)/plot["DY"].nBins),pred.GetBinContent(i))
        predErr.SetPointError(i, ((plot["DY"].firstBin-plot["DY"].lastBin)/plot["DY"].nBins)*0.5, ((plot["DY"].firstBin-plot["DY"].lastBin)/plot["DY"].nBins)*0.5, errFull, errFull)
        
        #up and down for ratio plot
        predUp.SetBinContent(i, pred.GetBinContent(i)+errSyst)
        predDown.SetBinContent(i, pred.GetBinContent(i)-errSyst)
        
    
    # Draw all kinds of things
    predErr.SetFillColor(myColors["MyBlue"])
    predErr.SetFillStyle(3001)
    
    data.UseCurrentStyle()
    
    pred.SetLineColor(ROOT.kBlack)
    pred.SetLineWidth(2)
    
    DYpred.SetLineColor(ROOT.kGreen+3)
    DYpred.SetFillColor(ROOT.kGreen+3)
    DYpred.SetFillStyle(3002)
    
    template.setPrimaryPlot(data, "PE")
    template.addSecondaryPlot(predErr, "02")
    template.addSecondaryPlot(pred, "hist")
    template.addSecondaryPlot(DYpred, "hist")
    
    template.hasRatio = True
    template.denominators = [pred]
    template.addRatioErrorByHist( "Syst. Error", predUp, predDown,color= myColors["MyBlue"],fillStyle=3001)           
    template.ratioLabel = "Data / Bgnd"
    
    template.dilepton = dil
    template.regionSize/=0.7
    template.draw()
    
    template.plotPad.cd()
    
    lines = getLines(0, data.GetMaximum()*0.6)
    for l in lines:
        l.Draw("same")
    
    leg = TLegend(0.6, 0.5, 0.9, 0.8)
    leg.SetFillColor(10)
    leg.SetTextFont(42)
    leg.SetLineColor(10)
    leg.SetShadowColor(0)
    leg.SetBorderSize(1)
    leg.AddEntry(data,"Data","pe")
    leg.AddEntry(pred,"Background Prediction","l")
    leg.AddEntry(DYpred,"Z+jets","f")
    leg.AddEntry(predErr,"Total Uncertainty","f")
    leg.Draw("same")
    
    template.ratioPad.cd()
    leg2 = TLegend(0.175, 0.78, 0.475, 0.9,"","brNDC")
    leg2.SetTextFont(42)
    leg2.SetFillColor(10)
    leg2.SetLineColor(10)
    leg2.SetShadowColor(0)
    leg2.SetBorderSize(1)
    leg2.AddEntry(predErr,"Syst. uncertainty", "f")   
    leg2.Draw("same")
        
    template.setFolderName("predictionMll")
    template.saveAs("spectrum_%s_%s"%(region, dil))
    template.clean()
    
def startPrediction(region):
    
    #initialize plots
    TT_full = dataMCConfig.dataMCConfig(plot="mllPlot",region="Signal"+region,runName="Run2015_25ns",plotData=True)
    DYScale = dataMCConfig.dataMCConfig(plot="mllPlotROutIn",region=getattr(regionsToUse.rOutIn, region.lower()).name,runName="Run2015_25ns",plotData=True)
    DY_full = dataMCConfig.dataMCConfig(plot="mllPlot",region=getattr(regionsToUse.rOutIn, region.lower()).name,runName="Run2015_25ns",plotData=True)
    
    plot = {"TT": TT_full.plot,
            "DYScale": DYScale.plot,
            "DY": DY_full.plot}

    #histograms for drell-yan prediction
    histDYEE, histDYMuMu, histDYOF = produceHistograms(plot["DY"])
    histDYSF = histDYEE.Clone()
    histDYSF.Add(histDYMuMu,1)
    
    #histograms for scaling of drell-yan spectra
    histDYScaleEE, histDYScaleMuMu, histDYScaleOF = produceHistograms(plot["DYScale"])
    histDYScaleSF = histDYScaleEE.Clone()
    histDYScaleSF.Add(histDYScaleMuMu,1)
                     
    #signal data + OF data for ttbar prediction
    dataEE, dataMuMu, histTT = produceHistograms(plot["TT"])
    dataSF = dataEE.Clone()
    dataSF.Add(dataMuMu,1)
    
    
    from corrections import rSFOF, rEEOF, rMMOF
    RSFOF = (getattr(rSFOF,region.lower()).val,getattr(rSFOF,region.lower()).err)
    REEOF = (getattr(rEEOF,region.lower()).val,getattr(rEEOF,region.lower()).err)
    RMMOF = (getattr(rMMOF,region.lower()).val,getattr(rMMOF,region.lower()).err)
    
    #run predictions on different dilepton combinations
    makePrediction(region, "SF",   plot, TT_full, RSFOF, dataSF,   histTT, histDYSF,   histDYScaleSF   )
    makePrediction(region, "EE",   plot, TT_full, REEOF, dataEE,   histTT, histDYEE,   histDYScaleEE   )
    makePrediction(region, "MuMu", plot, TT_full, RMMOF, dataMuMu, histTT, histDYMuMu, histDYScaleMuMu )
    
    
    
def main():
    #colors = createMyColors()   
    regions = ["Inclusive", "Forward", "Central"]
    #dataMCConfig.dataMCConfig.jzbType = "unCorrMet"
    for region in regions:
        startPrediction(region)
        
if __name__ == "__main__":
    main()
