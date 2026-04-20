import json
import subprocess
import os
import matplotlib.pylab as plt
import numpy as np

PATH = ""

ALLGATHER_TYPE = 0
ALLREDUCE_TYPE = 1
VERTEXSET_TYPE = 2

f = open(PATH + "experiments/test.json")
setup = json.load(f)
f.close()

experimentName = setup["experiment"]
nNodes = setup["nNodes"]
perNode = setup["perNode"]
mpilibs = setup["mpilib"]

# check for data locally
isLocal = True
for mpilib in mpilibs:
    isLocal = os.path.exists(PATH + "../data/{0}/{1}x{2}{3}.txt".format(experimentName,
                                                                        nNodes,
                                                                        perNode,
                                                                        mpilib))
    if not isLocal:
        break

if isLocal:
    ans = ""
    while ans not in ["y", "n", "Y", "N"]:
        ans = input("There is local data fetched. Do you want to use local data? (y/n)")
    if ans in ["n", "N"]:
        isLocal = False

if not isLocal:
    # retrieve data
    print("Processing experiment " + experimentName)
    print("Fetching data from cluster...")
    result = subprocess.run(
        "scp -r hydra:~/data/{0} {1}../data".format(experimentName, PATH),
        shell=True, capture_output=True, text=True
    )
    print("...done")

# parse data
graphs = setup[setup["plotGraphs"]]
xAxistype = setup["xAxis"]
xAxis = setup[xAxistype]


# sorting by graph
data = {g: {x: [] for x in xAxis} for g in graphs}

# data from different MPI libs are saved in different files 
if setup["plotGraphs"] == "mpilib":
    for g in graphs:
        f = open(PATH + "../data/{0}/{1}x{2}{3}.txt".format(experimentName,
                                                                        nNodes,
                                                                        perNode,
                                                                        g))
        line = f.readline()
        while line != "":
            chunks = line.split(",")
            mean = float (chunks[6].split(":")[1])
            if xAxistype == "setsize":
                xval = int (chunks[3].split(":")[1])            
            if xAxistype == "filling":
                xval = float (chunks[4].split(":")[1])
            data[g][xval].append(mean)
            line = f.readline()
        f.close()

# process and plot data
fig, ax = plt.subplots()
for g in graphs:
    # calculate error and mean for every x entry
    y = []
    yerr = []
    for x in xAxis:
        vals = np.array(data[g][x])
        y.append(np.mean(vals))
        yerr.append(np.std(vals))
    ax.errorbar(xAxis, y, yerr=yerr, ls="--", label=g, capsize=5)

ax.set_xscale("log")
ax.set_yscale("log")
ax.grid(True, which='major', linewidth=1.5)
ax.grid(True, which='minor', linewidth=0.5, linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()
