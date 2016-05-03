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

    args = parser.parse_args()
    
    if args.variables == []:
        from defs import varToUse
        args.variables = varToUse["JZB"].keys()
    
    ensurePathExists("shelves/")
    for var in args.variables:
        with open("shelves/%s.pkl"%(var), "w") as outputFile:
            pickle.dump(corrections, outputFile)
            outputFile.close()
        print "initialized shelves/%s.pkl"%(var)
        
if __name__ == "__main__":
    main()
