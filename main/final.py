import sys
sys.path.append("../")
from utils.final_check import *
from params.Funk_hist import *

if len(sys.argv) > 4:
	pathWithRedions = sys.argv[1]
	pathWithSource = sys.argv[2]
    pathToResults = sys.argv[3]
    pathToOut = sys.argv[4]
else:
	pathToRedions = "../../../Desktop/regions.txt"
    pathToSource = "../sirius_out/"
    pathToResults = "../../../Desktop/out/"
    pathToOut = "../../../Desktop/result.txt"


final_check(pathWithRedions,pathWithSource,pathToResults,pathToOut)

histo(pathToOut,)
