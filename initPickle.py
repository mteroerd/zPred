import pickle
import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

from helpers import ensurePathExists

def main():
    #corrections[Data][Correction]
    # soon corrections[Data][Direction][Correction]
    corrections = {}
    for plotData in [True, False]:
        corrections[plotData] = {}
        for direction in ["Central", "Forward","Inclusive"]:
            corrections[plotData][direction] = {}
            
            corrections[plotData][direction]["peak"] = {}
            corrections[plotData][direction]["pu"] = {}
            
            corrections[plotData][direction]["response"] = 0.0
            
            corrections[plotData][direction]["pu"][False] = 0.0
            corrections[plotData][direction]["pu"][True] = 0.0
            
            corrections[plotData][direction]["peak"][False] = 0.0
            corrections[plotData][direction]["peak"][True] = 0.0
    
    import argparse

    parser = argparse.ArgumentParser(description='Process some integers.')
        
    parser.add_argument("-v", "--vars", dest="variables", action="append", default=[],
                          help="JZB types to use")
    parser.add_argument("-c", "--corrs", dest="correctionMode", type=int, action="store", default=-1,
                          help="Which correction mode to use")
    args = parser.parse_args()
    
    
    if args.variables == []:
        from defs import varToUse
        args.variables = varToUse["JZB"].keys()
    
    ensurePathExists("shelves/")
    for var in args.variables:
        if args.correctionMode == -1:
            for i in range(0,4):
                with open("shelves/%s_%d.pkl"%(var,i), "w") as outputFile:
                    pickle.dump(corrections, outputFile)
                    outputFile.close()
                print "initialized shelves/%s_%d.pkl"%(var,i)
        else:
            with open("shelves/%s_%d.pkl"%(var,args.correctionMode), "w") as outputFile:
                pickle.dump(corrections, outputFile)
                outputFile.close()
            print "initialized shelves/%s_%d.pkl"%(var,args.correctionMode)
        
if __name__ == "__main__":
    main()
