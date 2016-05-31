import sys
import os.path
import os
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np
import argparse 
import time

from puCorr import doPUCorrection 
from peakCorr import doPeakCorrection
from responseCorr import doResponseCorrection

def doCorrs(var, args, extraArg):
    import dataMCConfig
    if not args.quiet: print "Calculating corrections for JZB type: %s"%(var)

    if not os.path.isfile("shelves/%s_%d.pkl"%(var,args.correctionMode)):
        print "Shelve file shelves/%s_%d.pkl does not exist yet, creating..."%(var,args.correctionMode)
        os.system("python initPickle.py --vars %s --corrs %d"%(var, args.correctionMode))
            
    if args.quiet or args.parallel:
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
    dataMCConfig.dataMCConfig.jzbType = var
    dataMCConfig.dataMCConfig.correctionMode = args.correctionMode
    
    if args.correctionMode > 2:
        if not args.quiet: print "%s: Response corrections..." %(var)
        for direction in ["Central", "Forward"]:
            if args.data:
                doResponseCorrection(True, direction, extraArg=extraArg)
            if args.mc:
                doResponseCorrection(False, direction, extraArg=extraArg)
        if not args.quiet: print "%s: Response corrections done." %(var)
    
    if args.correctionMode > 1:
        if not args.quiet: print "%s: Pile-up peak corrections..."%(var)
        for direction in ["Central", "Forward"]:
            if args.data:
                doPUCorrection(True, 2, direction, extraArg=extraArg)
                doPUCorrection(True, 3, direction, extraArg=extraArg)
            if args.mc:    
                doPUCorrection(False, 2, direction, extraArg=extraArg)
                doPUCorrection(False, 3, direction, extraArg=extraArg)
        if not args.quiet: print "%s: Pile-up peak corrections done."%(var)
    
    if args.correctionMode > 0:
        if not args.quiet: print "%s: Peak corrections..." %(var)
        for direction in ["Central", "Forward"]:
            if args.data:
                doPeakCorrection(True, 2, direction, extraArg=extraArg)
                doPeakCorrection(True, 3, direction, extraArg=extraArg)
            if args.mc:
                doPeakCorrection(False, 2, direction, extraArg=extraArg)
                doPeakCorrection(False, 3, direction, extraArg=extraArg)
        if not args.quiet: print "%s: Peak corrections done." %(var)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Process some integers.')
    
    parser.add_argument("-v", "--vars", dest="variables", action="append", default=[],
                          help="JZB types to use")
    parser.add_argument("-c", "--corrs", dest="correctionMode", type=int, action="store", default=-1,
                          help="Which correction mode to use")
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet", default=False,
                              help="shut up") 
    parser.add_argument("-p", "--parallel", action="store_true", dest="parallel", default=False,
                              help="do everything at once") 
    parser.add_argument("-d", "--data", action="store_true", dest="data", default=False,
                              help="calculate data corrections") 
    parser.add_argument("-m", "--mc", action="store_true", dest="mc", default=False,
                              help="calculate MC corrections") 
    
                              
    args = parser.parse_args()
    
    if args.correctionMode == -1:
        args.correctionMode = 3
    
    if (not args.data) and (not args.mc):
        args.data = True
        args.mc = True
        
    if args.variables == []:
        from defs import varToUse
        args.variables = varToUse["JZB"].keys()
    
    extraArg = "Q" if (args.quiet or args.parallel) else ""
     
    print "Starting"
    t1 = time.time()
    if not args.parallel:
        for var in args.variables:
            doCorrs(var, args, extraArg)
    else:
        import multiprocessing as mp
        processes = [mp.Process(target=doCorrs, args=(var, args, extraArg)) for var in args.variables]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
    t2 = time.time()    
    print "Finished."
    print "needed %f seconds"%(t2-t1)

if __name__ == "__main__":
    main()
