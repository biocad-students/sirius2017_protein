import sys
sys.path.append("../")
from utils.final_check import final_check
from params.Funk_hist import *

if len(sys.argv) > 1:
	pathToRegions = sys.argv[1] + "regions.txt"
	pathToSource = sys.argv[1]
else:
	pathToRegions = "../../../Desktop/sirius_out/regions.txt"
	pathToSource = "../../../Desktop/sirius_out/"

pathToRegions = "../../../Desktop/sirius_out/regions.txt"
pathToResults = "result/"
pathToOut = "rmsout.txt"


final_check(pathToRegions,pathToSource,pathToResults,pathToOut)
histo(pathToOut,0.1)
