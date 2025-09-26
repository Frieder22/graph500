import argparse
import sys
import ctypes

# Possible selections for modes and funcs
modes_selection = ["test", "perf"]
funcs_selection = ["EXACT_HALFING", "APROXIMATE_HALFING", "DENSE", "NAIVE"]
mpi_selection = {"openmpi": "spack load openmpi@4.1.6",
                  "mpich": "module load mpich/4.1.2-gcc-12.1.0-r7lq47v",
                  "intelmpi": "module load mpi/intelmpi"}

# define Parser stuff
parser = argparse.ArgumentParser(prog="createSlurmFile.py",
                                 description="Creates a slurm file that can be submitted on the hydra cluster",
                                 epilog="")

parser.add_argument("-m", "--mode", type=str, default="test",
                 help="Toggles, if a test or a performance measurement should be performed")
parser.add_argument("-n", "--nNodes", type=int, default=4,
                 help="Defines the number of nodes, that should be used")
parser.add_argument("-p", "--perNode", type=int, default=32,
                 help="Defines the number of processes per node")
parser.add_argument("-f", "--func", type=str, default="EXACT_HALFING",
                 help="Defines, which function should be tested or measured")
parser.add_argument("-mpi", "--mpi", type=str, default="openmpi",
                 help="Defines, which MPI library should be used")

args = parser.parse_args()

if args.mode not in modes_selection:
    print("Argparse error: The mode '" + args.mode + "' is not in " + str(modes_selection))
    sys.exit()

if args.func not in funcs_selection:
    print("Argparse error: The functionality (func) '" + args.func + "' is not in " + str(funcs_selection))
    sys.exit()

if not (args.nNodes > 0 and args.nNodes <= 36):
    print("Argparse error: The number of nodes (nNodes) '" + args.func + "' is not between 0 and 37. The hydra cluster has only 36 nodes available.")
    sys.exit()

if not (args.perNode > 0 and args.perNode <= 32):
    print("Argparse error: The number of processers per node (perNode) '" + args.func + "' is not between 0 and 33. The hydra cluster has only nodes with 32 cores each.")
    sys.exit()

if not args.mpi in mpi_selection:
    print("Argparse error: The MPI library (mpi) '" + args.mpi + "' is not in " + str(mpi_selection))
    sys.exit()   


mode = args.mode
func = args.func
nNodes = args.nNodes
perNode = args.perNode
mpi = args.mpi

# do some shananigans to find the correct enum values from "testStructure.h"
testnumber = {}
with open("../vertexSet/testStructure.h", "r") as f:
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

# create a slurm .job file
with open("run.job", "w") as f:
    f.write("#! /bin/bash\n")
    f.write("#SBATCH -p q_thesis\n")
    f.write("#SBATCH -N {}\n".format(nNodes))
    f.write("#SBATCH --ntasks-per-node={}\n".format(perNode))
    f.write("#SBATCH --cpu-freq=High\n")
    f.write("#SBATCH --time=3:00\n")
    if (mode == "perf"):
        f.write("srun ../build/{}.o {} {}\n".format(mode, testnumber[func], mpi))
    else:
        f.write("srun ../build/{}.o {} \n".format(mode, testnumber[func]))

# set correct mpi library in config file
with open("default.config", "w") as f:
    f.write("MPI_LOAD = {}\n".format(mpi_selection[mpi]))
    f.write("MPI_Variant = {}\n".format(mpi))