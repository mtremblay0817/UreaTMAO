import sys
import argparse
import MDAnalysis as mda
import numpy as np
from myTools.misc import readSimpleInput
from freud.locality import Voronoi
from time import perf_counter_ns
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

    required = ["top", "traj", "nFramesPerJob", "outStem", "waterName", "agent", "out"]

    # check for a bunch of possible errors
    errors = []
    for item in required:
        try:
            vars(args)[item]
        except KeyError:
            errors.append(f"ERROR: missing argument '{item}'")
        
    if errors:
        sys.exit("\n".join(errors))

    try:
        if jobIdx is None:
            args.start = 0
        else:
            args.out = f"{jobIdx}.results_temp"
            args.start = jobIdx * args.nFramesPerJob
        args.stop = args.start + args.nFramesPerJob
    except AttributeError:
        sys.exit("input file missing nFramesPerJob")

    frameDictionary = getDistribution(**vars(args))
    
    writeFrameData(frameDictionary, **vars(args))            

    print("Done <3")
    #pickAQuote("./quotes.txt")


def getDistribution(top, traj, waterName, agent, stop, start = 0,step = 1, **kwargs):
    """
    Finds the distribution of additions within the first and second shells in an MD simulation

    Arguments
    ---------
    top: str
        the topology file
    traj: str
        the trajectory file
    waterName: str
        the name of the water in the topology (SOL, HOH, etc.)
    agent: str
        the name of the addition (TMO, URE, etc.)
    start: int
        the first frame to analyze
    step: int
        the size of step between frames

    Returns
    -------
    agentCountFrames: dict[str : list[int]]
        the count of the agent in question and frames with that count
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

    protein = uni.select_atoms(f"not resname {waterName} NA CL Na+ Na Cl Cl- TMO URE UR", updating=True)
    proteinIdx = set(protein.ix)

    others = uni.select_atoms(f"resname {waterName} NA CL Na+ Na Cl Cl- TMO URE UR", updating=True)
    othersIdx = set(others.ix)
    residueMap = {atom.ix : atom.residue.ix for atom in others}

    agentCountFrames = {i : [] for i in range(100)}

    # nAtoms = {"TMO" : 14, "UR" : 8, "SOL" : 4, "CL" : 1}
    for ts in uni.trajectory[start:stop:step]:
        voronoiDiagram = Voronoi().compute(ts)
        neighborList = voronoiDiagram.nlist[:]

        firstAtoms = set()
        

        for i in range(len(neighborList)):
            if neighborList[i, 0] in proteinIdx and neighborList[i, 1] in othersIdx:
                firstAtoms.add(neighborList[i, 1])
        firstResidues = set()
        close = uni.select_atoms(f"byres (resname {waterName} NA CL Na Cl Na+ Cl- TMO URE UR and around 5.0 (not resname {waterName} NA CL Na Cl Na+ Cl- TMO URE UR))", periodic=True).ix
        for ix in np.intersect1d(close, list(firstAtoms)):
            firstResidues.add(residueMap[ix])
        firstShellResidues = uni.residues[list(firstResidues)].resnames
        agentCount = np.count_nonzero(firstShellResidues == agent)
        agentCountFrames[agentCount].append(ts.frame)

    return agentCountFrames

def writeFrameData(frameDictionary, out, **kwargs):
    """
    Writes the distribution of additions to an outfile

    Arguments
    ---------
    frameDictionary: dict[str : list[int]]
        the count of the agent in question and frames with that count
    framesOut: str
        the name of the outfile to write to

    Returns
    -------
    None
    """
    count_list = []
    with open(out, "w+") as outfile:
        for count in frameDictionary:
            if frameDictionary[count]:
                print(f"{count}:"+ ", ".join(list(map(str, frameDictionary[count]))), file=outfile)
                print(len(frameDictionary[count]), file = outfile)
                count_list += list(map(lambda x: count, frameDictionary[count]))
        print("Mean:", np.mean(count_list), file=outfile)
        print("StDev:", np.std(count_list), file=outfile)
        





if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("jobIdx", type=int, default=None, nargs="?")  # for graceParallel
    main(**vars(ap.parse_args()))
