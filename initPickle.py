import pickle
import sys
sys.path.append('cfg/')
from frameworkStructure import pathes
sys.path.append(pathes.basePath)

from helpers import ensurePathExists

def main():
    #corrections[Data][Correction]
    corrections = {True:{}, False:{}}
    corrections[True]["response"] = 0.0
    corrections[True]["pu"] = {}
    corrections[True]["pu"][False] = 0.0
    corrections[True]["pu"][True] = 0.0
    corrections[True]["peak"] = {}
    corrections[True]["peak"][False] = 0.0
    corrections[True]["peak"][True] = 0.0

    corrections[False]["response"] = 0.0
    corrections[False]["pu"] = {}
    corrections[False]["pu"][False] = 0.0
    corrections[False]["pu"][True] = 0.0
    corrections[False]["peak"] = {}
    corrections[False]["peak"][False] = 0.0
    corrections[False]["peak"][True] = 0.0

    import argparse

    parser = argparse.ArgumentParser(description='Process some integers.')
        
    parser.add_argument("-v", "--vars", dest="variables", action="append", default=[],
                          help="JZB types to use")
    parser.add_argument("-s", "--standard", dest="standard", action="store_true", default=False)
    parser.add_argument("-S", "--shift", dest="shift", action="store_true", default=False)
    args = parser.parse_args()
    
    if args.shift == False and args.standard == False:
        args.shift = True
        args.standard = True
    
    if args.variables == []:
        from defs import varToUse
        args.variables = varToUse["JZB"].keys()
    
    ensurePathExists("shelves/")
    for var in args.variables:
        if args.standard:
            with open("shelves/%s.pkl"%(var), "w") as outputFile:
                pickle.dump(corrections, outputFile)
                outputFile.close()
            print "initialized shelves/%s.pkl"%(var)
        if args.shift:
            with open("shelves/%s_onlyShift.pkl"%(var), "w") as outputFile:
                pickle.dump(corrections, outputFile)
                outputFile.close()
            print "initialized shelves/%s_onlyShift.pkl"%(var)
        
if __name__ == "__main__":
    main()
