import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)


import dataMCConfig

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import ratios
from defs import defineMyColors
from defs import myColors   

from setTDRStyle import setTDRStyle
from helpers import *   
colors = createMyColors() 
import math

################## SUMMARY OF CLASS plotTemplate ############################################################
## Constructors:
# * plotTemplate()
# * plotTemplate(mainConfig)
#       Applies some options from the mainConfig to the template:
#           axis titles from plot, private work, plot data, jzb type, integrated lumi from run range
#
## Methods:
# * setPrimaryPlot(plot, drawOption)
#       Plot (histogram, graph) to be drawn first with drawOption, defining the axes and boundaries, redrawn 
#       after secondary plots by default so it will be on top
#
# * addSecondaryPlot(plot, drawOption)
#       Adds plot (anything with Draw(drawOption) method in root) to list of secondary plots
#       In ratio graphs, primary plot / secondary plots is drawn, if denominators was not set
#       In efficiency graphs, secondary plots / primary plot is drawn
#
# * clearSecondaryPlots()
#       Resets list of secondary plots
#
# * addRatioErrorByHist(self,title, histUp, histDown, color, fillStyle)
#       See the RatioGraph class
#
# * addRatioErrorBySize(self,title, size, color, fillStyle, add)
#
#       See the RatioGraph class
# * addDenominator(histogram)
#       Add histogram to list of denominators for ratio plot by which mainPlot should be divided
#
# * clearDenominators()
#       Empty list of denominators
#
# * addNominator(histogram)
#       Add histogram to list of nominators for efficiency plot to be divided by mainPlot
#
# * clearNominators()
#       Empty list of nominators
#
# * addRatioPair(nominator, denominator, color)
#       Add a pair of histograms to be drawn in the ratioGraph. This way, secondaryPlots and (de)nominators will
#       be ignored for the ratio/efficiency and only ratioPairs will be drawn. Useful if there
#       is no real primary histogram to relate other histograms to, but multiple ratios are to be made or 
#       the ratios are not related to the plots in plotPad. Will also be used for efficiency plots.
#
# * clearRatioPairs()
#       Empty list of ratioPairs
#
# * draw()
#       Draw full canvas
#
# * drawLatexLabels()
#       If any of the labels were changed after calling draw(), this can be called to redraw all the labels
#
# * setFolderPath(folderName)
#       Set Name of folder in fig/ to store output
#
# * saveAs(fileName)
#       Print canvas to fileName in folder that was defined earlier, prints with all defined filetypes, so 
#       fileName does not contain file ending
#
# * clean()
#       Sets all objects in plotTemplate to None
#
#
## Members for options:
# -logX,logY,logZ(bools): 
#       Draw axis with logarithmic scale, false by default
# -changeScale(bool):
#       True by default, defines if the minimum/maximum of primaryPlot should be adjusted by 
#       maximumScale or minimum/maximum. Should be turned to False if primaryPlot is not a Histogramm
# -maximumScale(float):
#       Maximum of plot scaled by maximum value of primary plot
# -minimum, maximum(floats):
#       Overwrites maximum scale, set minimum and maximum y-value(z-value for 2D plots) of primary plot
# -labelX,labelY,labelZ(strings): 
#       Set axis title of primary plot, None (default) will not change titles of primary plot
# -marginTop, marginBottom, marginLeft, marginRight(floats): 
#       Pad margins
# -jzbType, useJZBPath(string, bool): 
#       Type of used jzb, should it be used in the path when calling saveAs (True by default)
# -personalWork, preliminary, forTWIKI, plotData(bools): 
#       For text below CMS logo. Default: True, False, False, False
# -hasRatio, hasEfficiency(bools): 
#       Draw ratio/efficiency graph, if both are true, only ratio graph is drawn. Default: False
# -ratioLabel, efficiencyLabel(strings): 
#       Y-axis label of ratio/efficiency graph
# -ratioMin, ratioMax(floats):
#       Min/Max values in ratio graph
# -denominators(list of histograms):
#       List of plots to divide primary plot in ratio graph. If empty, all secondary plots will be used
# -nominators(list of histograms):
#       List of plots to be divided by primary plot in efficiency graph. If empty, all secondary plots will 
#       be used
# -ratioPairs (list of tuples of histogram, histogram, color):
#       If not equal to [], will ignore denominators and secondaryPlots and add ratios of the given pairs to
#       the ratioPlot. Will also be used for efficiency plots.
# -redrawPrimary(bool):
#       Should primary plot be redrawn to be on top of other plots. On by default, not recommended for 2D plots
# -dilepton (string):
#       Used dilepton combination in plot. "SF", "EE", "EMu", "MuMu" will result in predefined dilepton strings.
#       If dilepton is set to a different string, this will be used directly.
# -fileTypes (list of strings):
#       List of file endings to print the canvas into. Default: ["pdf"]
#  ---
# Following members are for latex labels of cms, cmsExtra, region, cuts
#  -#PosX, #PosY (floats): 
#       Position of label in NDC. lumiPosX is always used with -(marginRight-0.05) as an offset, so it
#       will be aligned to the edge of the plot if a right margin is introduced
#  -cmsExtraSimPosY (float):
#       Position of cmsExtra label if using simulated data in private work
#  -#Size (float):
#       Text size
#  -#Font, #Align (ints):
#       Text font and align
#  -#Do (bool):
#       Should label be printed (not defined for cuts)
#  -cutsText(string): 
#       Text to describe cuts; if not set, the label will not be drawn
#  -lumiInt, lumiSqrtS (floats):
#       Integrated luminosity [fb^-1] and sqrt(s) [TeV]; if lumiInt is not set, no luminosity will be printed

def countNumbersUp():
    countNumbersUp.counter += 1
    return countNumbersUp.counter
countNumbersUp.counter = 0

class plotTemplate:    
    useJZBPath = True
    pathName = "fig/"
    folderName = ""
    fileTypes = ["pdf"]
    dilepton = None
    
    logX = False
    logY = False
    logZ = False
    redrawPrimary = True
    maximumScale = 1.8
    maximum = None
    minimum = None
    changeScale = True
    
    marginLeft = 0.15
    marginRight = 0.05
    marginTop = 0.05
    marginBottom = 0.14
    
    ##############
    latexCMS = None
    cmsPosX = 0.19 
    cmsPosY = 0.89
    cmsSize = 0.06
    cmsFont = 61
    cmsAlign = 11
    cmsColor = ROOT.kBlack
    cmsDo = True
    ##############
    latexCMSExtra = None
    cmsExtraPosX = 0.19
    cmsExtraPosY = 0.85
    cmsExtraSimPosY = 0.82
    cmsExtraSize = 0.045
    cmsExtraFont = 52
    cmsExtraAlign = 11
    cmsExtraColor = ROOT.kBlack
    cmsExtraDo = True
    ##############
    latexLumi = None
    lumiPosX = 0.95 
    lumiPosY = 0.96
    lumiSize = 0.04
    lumiFont = 42
    lumiAlign = 31
    lumiSqrtS = 13
    lumiInt = None
    lumiColor = ROOT.kBlack
    lumiDo = True
    ##############
    latexRegion = None
    latexRegion2 = None
    regionPosX = 0.92
    regionPosY = 0.89
    regionSize = 0.03
    regionFont = 42
    regionAlign = 32
    regionColor = ROOT.kBlack
    regionName = ""
    regionDo = True
    ##############
    latexCuts = None
    cutsPosX = 0.92
    cutsPosY = 0.81
    cutsSize = 0.03
    cutsFont = 42
    cutsAlign = 32
    cutsText = None
    cutsColor = ROOT.kBlack
    ##############
    
    labelX = None
    labelY = None
    labelZ = None 
    
    ratioLabel = "ratio"    
    efficiencyLabel = "efficiency"
    ratioMin = 0
    ratioMax = 2
    ratioErrsSize = []
    ratioErrsHist = []
    denominators = []
    nominators = []
    ratioPairs = []
    hasEfficiency = False
    hasRatio = False
    
    jzbType = "type-IMet"
    
    forTWIKI = False
    preliminary = False
    personalWork = True
    plotData = False
    
    def __init__(self,mainConfig=None):        
        self.canvas = None
        self.plotPad = None
        self.ratioPad = None
        
        self.primaryPlot = None
        self.secondaryPlots = []
        
        self.ratioErrsSize = []
        self.ratioErrsHist = []
        
        self.ratioPairs = []
        self.nominators = []
        self.denominators = []
        setTDRStyle() 
        
        self.latexCMS      = None
        self.latexCMSExtra = None
        self.latexCuts     = None
        self.latexLumi     = None
        self.latexRegion   = None
        
        if mainConfig != None:
            self.regionName = mainConfig.plot.label2
            if hasattr(mainConfig, "jzbType"):
                self.jzbType = mainConfig.jzbType
            self.personalWork = mainConfig.personalWork
            self.preliminary = mainConfig.preliminary
            self.forTWIKI = mainConfig.forTWIKI
            self.plotData = mainConfig.plotData
            self.hasRatio = mainConfig.plotRatio
            self.lumiInt = mainConfig.runRange.printval
            self.labelX = mainConfig.plot.xaxis
            self.labelY = mainConfig.plot.yaxis
    
    def addRatioPair(self, nominator, denominator, color=ROOT.kBlack):
        self.ratioPairs.append((nominator, denominator, color))
        
    def clearRatioPairs(self):
        self.ratioPairs = []
        
    def addNominator(self, nominator):
        self.nominators.append(nominator)
        
    def clearNominators(self):
        self.nominators = []
        
    def addDenominator(self, denominator):
        self.denominators.append(denominator)
    
    def clearDenominators(self):
        self.denominators = []

    def addRatioErrorBySize(self,title, size, color, fillStyle, add):
        self.ratioErrsSize.append((title,size,color,fillStyle,add))

    def addRatioErrorByHist(self,title, histUp, histDown, color, fillStyle):
        self.ratioErrsHist.append((title,histUp,histDown,color,fillStyle))
        
    def setFolderName(self,name):
        self.folderName = name
        if self.useJZBPath:
            self.pathName = "fig/%s/%s/"%(self.jzbType, name)
        else:
            self.pathName = "fig/%s/"%(name)
    
    def saveAs(self,fileName):
        ensurePathExists(self.pathName)
        for typ in self.fileTypes:
            self.canvas.Print(self.pathName+fileName+"."+typ)
    
    def setPrimaryPlot(self,hist, drawOption):
        self.primaryPlot = (hist, drawOption)
    
    def addSecondaryPlot(self,hist, drawOption=""):
        self.secondaryPlots.append((hist,drawOption))
    
    def clearSecondaryPlots(self):
        self.secondaryPlots = []
        
    def drawLatexLabels(self):
        #CMS Text
        if self.cmsDo:
            if self.latexCMS == None:
                self.latexCMS = ROOT.TLatex()
            self.latexCMS.SetTextFont(self.cmsFont)
            self.latexCMS.SetTextSize(self.cmsSize)
            self.latexCMS.SetTextAlign(self.cmsAlign)
            self.latexCMS.SetNDC(True)
            self.latexCMS.SetTextColor(self.cmsColor)
            self.latexCMS.DrawLatex(self.cmsPosX,self.cmsPosY,"CMS")
        else:
            self.latexCMS = None
        #Sub CMS Text
        if self.cmsExtraDo:
            if self.latexCMSExtra == None:
                self.latexCMSExtra = ROOT.TLatex()
            self.latexCMSExtra.SetTextFont(self.cmsExtraFont)
            self.latexCMSExtra.SetTextSize(self.cmsExtraSize)
            self.latexCMSExtra.SetTextAlign(self.cmsExtraAlign)
            self.latexCMSExtra.SetNDC(True) 
            self.latexCMSExtra.SetTextColor(self.cmsExtraColor)
            yLabelPos = self.cmsExtraPosY
            cmsExtra = ""
            if self.personalWork:
                cmsExtra = "Private Work"
                if not self.plotData:
                    cmsExtra = "#splitline{Private Work}{Simulation}"
                    yLabelPos = self.cmsExtraSimPosY
            elif not self.plotData:
                cmsExtra = "Simulation" 
            elif self.preliminary:
                cmsExtra = "Preliminary"
            elif self.forTWIKI:
                cmsExtra = "Unpublished"        
            self.latexCMSExtra.DrawLatex(self.cmsExtraPosX,yLabelPos,"%s"%(cmsExtra))
        else:
            self.latexCMSExtra = None
        #Lumi and sqrt(s) Text
        if self.lumiDo:
            if self.latexLumi == None:
                self.latexLumi = ROOT.TLatex()
            self.latexLumi.SetTextFont(self.lumiFont)
            self.latexLumi.SetTextAlign(self.lumiAlign)
            self.latexLumi.SetTextSize(self.lumiSize)
            self.latexLumi.SetNDC(True)  
            self.latexLumi.SetTextColor(self.lumiColor)                      
            if self.lumiInt != None:
                self.latexLumi.DrawLatex(self.lumiPosX-(self.marginRight-0.05), self.lumiPosY, "%s fb^{-1} (%s TeV)"%(self.lumiInt,self.lumiSqrtS))   
            else:
                self.latexLumi.DrawLatex(self.lumiPosX-(self.marginRight-0.05), self.lumiPosY, "%s TeV"%(self.lumiSqrtS))
        else:
            self.latexLumi = None
        #Region identifier
        if self.regionDo:
            if self.latexRegion == None:
                self.latexRegion = ROOT.TLatex()
            self.latexRegion.SetTextAlign(self.regionAlign)
            self.latexRegion.SetTextSize(self.regionSize)
            self.latexRegion.SetNDC(True)
            self.latexRegion.SetTextFont(self.regionFont)
            self.latexRegion.SetTextColor(self.regionColor)
            if self.dilepton != None:
                if self.dilepton == "SF":
                    dileptonLabel = "ee+#mu#mu"
                elif self.dilepton == "EE":
                    dileptonLabel = "ee"
                elif self.dilepton == "MuMu":
                    dileptonLabel = "#mu#mu"
                elif self.dilepton == "EMu":
                    dileptonLabel = "e#mu"
                else:
                    dileptonLabel = self.dilepton
            else:
                dileptonLabel = ""
            
            if dileptonLabel != "":
                if self.regionName != "":
                    self.latexRegion2 = self.latexRegion.Clone()
                    self.latexRegion.DrawLatex(self.regionPosX,self.regionPosY+0.5*self.regionSize,self.regionName) 
                    self.latexRegion2.DrawLatex(self.regionPosX,self.regionPosY-0.5*self.regionSize,dileptonLabel)
                else:
                    self.latexRegion.DrawLatex(self.regionPosX,self.regionPosY,dileptonLabel) 
            else:
                self.latexRegion.DrawLatex(self.regionPosX,self.regionPosY,self.regionName) 
        else:
            self.latexRegion = None
        #Cuts
        if self.cutsText != None:
            if self.latexCuts == None:
                self.latexCuts = ROOT.TLatex()
            self.latexCuts.SetTextFont(self.cutsFont)
            self.latexCuts.SetTextAlign(self.cutsAlign)
            self.latexCuts.SetTextSize(self.cutsSize)
            self.latexCuts.SetNDC(True)   
            self.latexCuts.SetTextColor(self.cutsColor)    
            self.latexCuts.DrawLatex(self.cutsPosX, self.cutsPosY, self.cutsText)
        else:
            self.latexCuts = None
            
    def clean(self):       
        self.latexCMS       = None
        self.latexCMSExtra  = None
        self.latexCuts      = None
        self.latexLumi      = None
        self.latexRegion    = None
        self.primaryPlot    = None
        self.secondaryPlots = []
        self.ratioGraphs    = []
        self.denominators   = []
        self.nominators     = []
        self.ratioPairs     = []
        self.pathName       = "fig/"
        self.folderName     = ""
        self.canvas         = None
        self.plotPad        = None
        self.ratioPad       = None
    
    def draw(self):
        self.canvas = ROOT.TCanvas("hCanvas%d"%(countNumbersUp()), "", 800,800)
        if self.hasRatio or self.hasEfficiency:
            self.plotPad = ROOT.TPad("plotPad","plotPad",0,0.3,1,1)
        else:
            self.plotPad = ROOT.TPad("plotPad","plotPad",0,0,1,1)
        self.plotPad.UseCurrentStyle()
        self.plotPad.Draw()  
        
        if self.hasRatio or self.hasEfficiency:
            self.ratioPad = ROOT.TPad("ratioPad","ratioPad",0,0,1,0.3)
            self.ratioPad.UseCurrentStyle()
            self.ratioPad.Draw()
         
        if self.hasRatio or self.hasEfficiency:
            self.plotPad.SetTopMargin    (self.marginTop)
            self.plotPad.SetLeftMargin   (self.marginLeft)
            self.plotPad.SetRightMargin  (self.marginRight)
            self.ratioPad.SetBottomMargin(self.marginBottom)
            self.ratioPad.SetLeftMargin  (self.marginLeft)
            self.ratioPad.SetRightMargin (self.marginRight)
        else:
            self.plotPad.SetTopMargin   (self.marginTop)
            self.plotPad.SetLeftMargin  (self.marginLeft)
            self.plotPad.SetRightMargin (self.marginRight)
            self.plotPad.SetBottomMargin(self.marginBottom)
            
        self.plotPad.cd()  
        
        if self.logX:
            self.plotPad.SetLogx()
            if self.hasRatio or self.hasEfficiency:
                self.ratioPad.SetLogx()
        if self.logY:
            self.plotPad.SetLogy()
        if self.logZ:
            self.plotPad.SetLogz()
          
        if self.changeScale:
            if self.maximum == None:
                self.primaryPlot[0].SetMaximum(self.primaryPlot[0].GetMaximum()*self.maximumScale)
            else:
                self.primaryPlot[0].SetMaximum(self.maximum)
            if self.minimum != None:
                self.primaryPlot[0].SetMinimum(self.minimum)
        
        if self.labelX != None:
            self.primaryPlot[0].GetXaxis().SetTitle(self.labelX)
        if self.labelY != None:
            self.primaryPlot[0].GetYaxis().SetTitle(self.labelY)
        if self.labelZ != None:
            self.primaryPlot[0].GetZaxis().SetTitle(self.labelZ)
            
        self.primaryPlot[0].Draw(self.primaryPlot[1])
        
        for plot, drawStyle in self.secondaryPlots:
            plot.Draw(drawStyle+"same")
         
        if self.redrawPrimary:
            self.primaryPlot[0].Draw(self.primaryPlot[1]+"same")
         
        
        self.drawLatexLabels()
        self.plotPad.RedrawAxis()
        
        self.drawRatioPlots()
        
    def drawRatioPlots(self):
        if self.hasRatio or self.hasEfficiency:
            self.ratioPad.cd()
            if self.hasRatio:
                self.ratioGraphs = []
                if self.ratioPairs == []:
                    if self.denominators == []:
                        self.denominators = [plot[0] for plot in self.secondaryPlots]
                    for denominator in self.denominators:
                        self.ratioGraphs.append(ratios.RatioGraph(self.primaryPlot[0],denominator, xMin=self.primaryPlot[0].GetXaxis().GetBinLowEdge(1), xMax=self.primaryPlot[0].GetXaxis().GetBinUpEdge(self.primaryPlot[0].GetNbinsX()),title=self.ratioLabel,yMin=self.ratioMin,yMax=self.ratioMax,ndivisions=10,color=denominator.GetMarkerColor(),  adaptiveBinning=1000))
                    for err in self.ratioErrsSize:
                        self.ratioGraphs[0].addErrorBySize(err[0],err[1],err[2],err[3],err[4])
                    for err in self.ratioErrsHist:
                        self.ratioGraphs[0].addErrorByHistograms(err[0],err[1],err[2],err[3],err[4])
                else:
                    for nominator, denominator, color in self.ratioPairs:
                        self.ratioGraphs.append(ratios.RatioGraph(nominator, denominator, xMin=self.primaryPlot[0].GetXaxis().GetBinLowEdge(1), xMax=self.primaryPlot[0].GetXaxis().GetBinUpEdge(self.primaryPlot[0].GetNbinsX()),title=self.ratioLabel,yMin=self.ratioMin,yMax=self.ratioMax,ndivisions=10,color=color,  adaptiveBinning=1000 ))
                for number,graph in enumerate(self.ratioGraphs):
                    if number == 0:
                        graph.draw(ROOT.gPad,True,False,True,chi2Pos=0.8)
                    else:
                        graph.draw(ROOT.gPad,False,False,True,chi2Pos=0.8)
            elif self.hasEfficiency:
                self.ratioGraphs = []
                if self.nominators == []:
                    self.nominators = [plot[0] for plot in self.secondaryPlots]
                self.hAxis = ROOT.TH2F("hAxis%d"%(countNumbersUp()), "", 20, self.primaryPlot[0].GetXaxis().GetBinLowEdge(1), self.primaryPlot[0].GetXaxis().GetBinUpEdge(self.primaryPlot[0].GetNbinsX()), 10, 0, 1.2)    
                self.hAxis.GetYaxis().SetNdivisions(408)
                self.hAxis.GetYaxis().SetTitleOffset(0.4)
                self.hAxis.GetYaxis().SetTitleSize(0.15)
                self.hAxis.GetXaxis().SetLabelSize(0.0)
                self.hAxis.GetYaxis().SetLabelSize(0.15)
                self.hAxis.GetYaxis().SetTitle(self.efficiencyLabel)
                self.hAxis.Draw("AXIS")
                
                if self.ratioPairs == []:
                    for nominator in self.nominators:
                        tmp = ROOT.TGraphAsymmErrors(nominator,self.primaryPlot[0], "cp")
                        tmp.SetMarkerColor(nominator.GetMarkerColor())
                        tmp.SetLineColor(nominator.GetLineColor())
                        self.ratioGraphs.append(tmp)
                        self.ratioGraphs[len(self.ratioGraphs)-1].Draw("same P")
                else:
                    for nominator, denominator, color in ratioPairs:
                        tmp = ROOT.TGraphAsymmErrors(nominator,denominator, "cp")
                        tmp.SetMarkerColor(color)
                        tmp.SetLineColor(color)
                        self.ratioGraphs.append(tmp)
                        self.ratioGraphs[len(self.ratioGraphs)-1].Draw("same P")


class plotTemplate2D(plotTemplate):
    marginRight = 0.2
    regionPosX = 0.78
    cutsPosX = 0.78
    redrawPrimary = False
    maximumScale = 1.1
    
    def __init__(self, mainConfig=None):
        plotTemplate.__init__(self, mainConfig)
        if mainConfig is not None:
            if hasattr(mainConfig, "plot2"):
                if mainConfig.plot2 is not None:
                    self.labelY = mainConfig.plot2.xaxis
