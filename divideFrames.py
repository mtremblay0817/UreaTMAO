import sys
import argparse
import MDAnalysis as mda
import numpy as np
from myTools.misc import readSimpleInput
from freud.locality import Voronoi
from tqdm import tqdm
# from pickAQuote import pickAQuote


def main(*, infile, jobIdx=0):
    """
    Calculate distribution of additives in first hydration shells of proteins

    Arguments
    ---------
    infile: str
        input file name, must be in the correct format for readInputFile()
    jobIdx: int (default = 0)
        job number, 0 - (nJobs-1), (for parallel calculations)

    """
    print("\nProcessing input...")
    args = readSimpleInput(infile)

    required = ["top", "traj", "data", "outStem", "splits"]

    # check for a bunch of possible errors
    errors = []
    for item in required:
        try:
            vars(args)[item]
        except KeyError:
            errors.append(f"ERROR: missing argument '{item}'")
        
    if errors:
        sys.exit("\n".join(errors))

    # try:
    #     if jobIdx is None:
    #         args.start = 0
    #     else:
    #         args.out = f"{jobIdx}.results_temp"
    #         args.start = jobIdx * args.nFramesPerJob
    #     args.stop = args.start + args.nFramesPerJob
    # except AttributeError:
    #     sys.exit("input file missing nFramesPerJob")

    divideTraj(**vars(args))
    print("Done <3")
    #pickAQuote("./quotes.txt")


def divideTraj(top, traj, data, outStem, splits, **kwargs):
    """
    Finds the distribution of additions within the first and second shells in an MD simulation

    Arguments
    ---------
    top: str
        the topology file
    traj: str
        the trajectory file
    data: str
        the file containing the frame data
    outStem: str
        the name stem of the trajectory
    splits: str
        the division of counts into separate trajectories
        format - start1;stop1,start2;stop2....
        Can use ; with empty side like slices use :

    Returns
    -------
    None
    """

    # kill warning about unknown atom masses, useless
    mda.warnings.simplefilter("ignore")
    if top.endswith(".top"):
        try:
            print("reading GROMACS topology...")
            uni = mda.Universe(top, traj, topology_format="ITP")
        except Exception:
            print("GROMACS topology: make sure you have #included full paths to force field in topology!")
            raise RuntimeError(".top files are read as GROMACS topologies, change topology extension")
    else:
        print("reading structure/topology...")
        uni = mda.Universe(top, traj)

    all_atoms = uni.select_atoms("all")
    print("Parsing frame data", flush=True)
    datafile = open(data)
    frames = {}
    for line in datafile.readlines():
        n, f = line.split(":")
        frames[int(n)] = list(map(int, f.split(",")))
    datafile.close()
    print("Parsing trajectory splits", flush=True)
    if splits[0][0] == ";":
        splits[0] = str(min(frames.keys()))+splits[0]
    if splits[-1][-1] == ";":
        splits[-1] += str(max(frames.keys()))
    print("Combining splits", flush=True)
    trajSplits = {}

    for section in splits:
        start,stop = section.split(";")
        splitName = f"{start}_{stop}"
        trajSplits[splitName] = []
        for i in range(int(start), int(stop) + 1):
            trajSplits[splitName] += frames[i]

    print("Writing trajectories", flush=True)
    for splitName in tqdm(trajSplits):
        all_atoms.write(outStem[:-4]+splitName+".dcd", frames=trajSplits[splitName])


if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("jobIdx", type=int, default=None, nargs="?")  # for graceParallel
    main(**vars(ap.parse_args()))
