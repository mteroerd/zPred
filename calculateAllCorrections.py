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
from Response import plotResponse

def doCorrs(var, args, extraArg):
    import dataMCConfig
    if not args.quiet: print "Calculating corrections for JZB type: %s"%(var)
    if not os.path.isfile("shelves/%s.pkl"%(var)):
        print "Shelve file shelves/%s.pkl does not exist yet, creating..."%(var)
        os.system("python initPickle.py --vars %s"%(var))
    if args.quiet or args.parallel:
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
    dataMCConfig.dataMCConfig.jzbType = var
    if not args.quiet: print "%s: Response corrections..." %(var)
    if args.data:
        plotResponse("SF",True, extraArg=extraArg)
    if args.mc:
        plotResponse("SF",False, extraArg=extraArg)
    if not args.quiet: print "%s: Response corrections done." %(var)
    
    if not args.quiet: print "%s: Pile-up peak corrections..."%(var)
    if args.data:
        doPUCorrection("SF",True,2, extraArg=extraArg)
        doPUCorrection("SF",True,3, extraArg=extraArg)
    if args.mc:    
        doPUCorrection("SF",False,2, extraArg=extraArg)
        doPUCorrection("SF",False,3, extraArg=extraArg)
    if not args.quiet: print "%s: Pile-up peak corrections done."%(var)
    
    if not args.quiet: print "%s: Peak corrections..." %(var)
    if args.data:
        doPeakCorrection("SF", True,2, extraArg=extraArg)
        doPeakCorrection("SF", True,3, extraArg=extraArg)
    if args.mc:
        doPeakCorrection("SF", False,2, extraArg=extraArg)
        doPeakCorrection("SF", False,3, extraArg=extraArg)
    if not args.quiet: print "%s: Peak corrections done." %(var)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Process some integers.')
    
    parser.add_argument("-v", "--vars", dest="variables", action="append", default=[],
                          help="JZB types to use")
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet", default=False,
                              help="shut up") 
    parser.add_argument("-p", "--parallel", action="store_true", dest="parallel", default=False,
                              help="do everything at once") 
    parser.add_argument("-d", "--data", action="store_true", dest="data", default=False,
                              help="calculate data corrections") 
    parser.add_argument("-m", "--mc", action="store_true", dest="mc", default=False,
                              help="calculate MC corrections") 
    
                              
    args = parser.parse_args()
    
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
