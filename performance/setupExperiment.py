import numpy as np
import json
import random

#PATH = "performance/"     # used for using vsCODE debugger
PATH = ""

NEXPERIMENTRUNS = 30

MPI_LOAD = {"openMPI": "spack load openmpi@4.1.6",
           "intelMPI": "module load mpi/intelmpi",
           "MPICH": "module unload mpi/intelmpi"
           }

MPI_UNLOAD = {"openMPI": "spack unload openmpi@4.1.6 ",
           "intelMPI": "module unload mpi/intelmpi",
           "MPICH": "module unload mpi/intelmpi"
           }

MPI_LAUNCHER = {"openMPI": "srun",
           "intelMPI": "srun",
           "MPICH": "mpirun"
           }

def createRunsJSON(data, tn):
    # change all listable keys into lists
    for key in ["funcType", "vsFunc", "setSize", "filling"]:
        if type(data[key]) != list:
            data[key] = [data[key]]

    # find all run setups
    runList = []
    for funcType in data["funcType"]:
        for setSize in data["setSize"]:
            for filling in data["filling"]:
                if funcType == "Vertexset":
                    for vsFunc in data["vsFunc"]:
                        runList.append((funcType, vsFunc, setSize, filling))
                else:
                    runList.append((funcType, "ITERATOR", setSize, filling)) # using ITERATOR as dummy

    # shuffle runs
    random.seed(42)
    random.shuffle(runList)
    # create run.json file
    outJSON = {
        "runs": [
            {"experiment": data["experiment"],
            "funcType": funcType,
            "vsFunc": tn[vsFunc],
            "setSize": setSize,
            "filling": filling,
            "nRepetitions": data["nRepetitions"]}
            for funcType, vsFunc, setSize, filling in runList
        ]
    }
    f = open(PATH + "runs.json", "w")
    json.dump(outJSON, f, indent=4)
    f.close()

def createSlurmJob(data):
    # prepare and get setup values
    nNodes = data["nNodes"]
    perNode = data["perNode"]
    if type(data["mpilib"]) != list:
        data["mpilib"] = [data["mpilib"]]

    f = open("run.job", "w")
    # setup on cluster
    f.write("#! /bin/bash\n")
    f.write("#SBATCH -p q_thesis\n")
    f.write("#SBATCH -N {}\n".format(nNodes))
    f.write("#SBATCH --ntasks-per-node={}\n".format(perNode))
    f.write("#SBATCH --cpu-freq=High\n")
    f.write("#SBATCH --job-name=performanceRun\n")
    f.write("#SBATCH --time=1:00:00\n\n")

    #create folder for saved data
    f.write("mkdir ../data/{0}\n\n".format(data["experiment"]))
    # procedure for different MPI libs
    for mpilib in data["mpilib"]:
        # delete old measurement data, if there is some present
        f.write("rm ../data/{0}/{1}x{2}{3}.txt\n".format(data["experiment"],
                                                data["nNodes"],
                                                data["perNode"],
                                                mpilib))
        
        f.write("echo MPIlib: {0}\n".format(mpilib))
        f.write(MPI_LOAD[mpilib] +  "\n")
        f.write("make perf\n")
        f.write("for i in {{1..{0}}}; do\n".format(NEXPERIMENTRUNS))
        f.write("    {0} ../build/executeExperimentRuns.o runs.json >> ../data/{1}/{2}x{3}{4}.txt\n".format(MPI_LAUNCHER[mpilib],
                                                                                                     data["experiment"],
                                                                                                     data["nNodes"],
                                                                                                     data["perNode"],
                                                                                                     mpilib))
        f.write("echo iteration $i\n")
        f.write("done\n")
        f.write(MPI_UNLOAD[mpilib] + "\n\n")

    f.write("echo finished\n")

    f.close()

def findTestNumber():
    # do some shananigans to find the correct enum values from "testStructure.h"
    testnumber = {}
    with open(PATH + "../vertexSet/testStructure.h", "r") as f:
        line = f.readline()
        while line != "enum Testcase{\n":
            line = f.readline()
        count = 0
        line = f.readline()
        while line != "};\n":
            word = line.split(",")[0].replace(" ", "").replace("\n","")
            testnumber[word] = count
            count += 1
            line = f.readline()
    return testnumber


# load setup data
f = open(PATH + "experiments/test.json")
toExecuteData = json.load(f)
f.close()

# find correct IDs for Vertexsetfunctions
tn = findTestNumber()

createRunsJSON(toExecuteData, tn)
createSlurmJob(toExecuteData)
