import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from myTools.misc import readSimpleInput
from pickAQuote import pickAQuote


def main(infile):
    print("\nProcessing input...")
    args = readSimpleInput(infile)

    required = ["top", "traj", "restypes", "frames", "out", "plot"]

    # check for a bunch of possible errors
    errors = []
    for item in required:
        try:
            vars(args)[item]
        except AttributeError:
            errors.append(f"ERROR: missing argument '{item}'")


    if isinstance(args.restypes, str):
        args.restypes = [args.restypes.upper()]
    else:
        for i in range(len(args.restypes)):
            args.restypes[i] = args.restypes[i].upper()

    if args.longitudinal is True and len(args.restypes) > 1:
        errors.append("ERROR: Only one restype is allowed for longitudinal")
    elif args.plot.lower()=="mean" and args.longitudinal is False:
        errors.append("ERROR: Longitudinal is mandatory for mean curves")
        
    if errors:
        sys.exit("\n".join(errors))

    output = getDistribution(**vars(args))

    with open(args.out, "w+") as o:
        for add, dist in output.items():
            o.write(add+":"+", ".join(list(map(str, dist)))+"\n")

    if args.plot.lower()=="hist":
        print("Plotting...")
        plot_hist(output, **vars(args))
    elif args.plot.lower()=="mean":
        print("Plotting...")
        plot_mean(output, **vars(args))

    pickAQuote("./quotes.txt")


def getDistribution(top, traj, restypes, frames, frameStart = 0, frameStep = 1, longitudinal = False, **kwargs):
    """
    Finds the distribution of additions to an MD simulation

    Arguments
    ---------
    top: str
        the topology file
    traj: str
        the trajectory file
    restypes: list[str]
        the residues being examined
    frames: int or str
        total number of frames to analyze or "all"
    frameStart: str
        the first frame to analyze
    frameStep: str
        the size of step between frames
    longitudinal: bool
        True if examining one restype over multiple frames, otherwise accumulate over all frames

    returns
    -------
    addDistances: dict[str, list]
        distances of residues in question from nearest protein atom

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
    
    if type(frames) == int:
        frameStop = frameStart + frames*frameStep
    elif frames == "all":
        frameStop = len(uni.trajectory)
    else:
        sys.exit("Invalid frame selection")


    prot = uni.select_atoms('not resname TMO UR SOL CL', updating=True)

    # nAtoms = {"TMO" : 14, "UR" : 8, "SOL" : 4, "CL" : 1}

    if longitudinal:
        res = restypes[0]
        addDistances = {}
        sel = uni.select_atoms(f"resname {res}", updating=True).split('residue')
        
        for ts in uni.trajectory[frameStart:frameStop:frameStep]:
            print(ts.frame, flush=(((ts.frame-frameStart)/frameStep)%10==0))
            addDistances[f"{res} Frame {ts.frame}"] = []
            for residue in sel:
                dist_array = distances.distance_array(residue.positions, prot.positions,box=uni.dimensions)
                addDistances[f"{res} Frame {ts.frame}"].append(np.min(dist_array))
    else:
        addDistances = {res : [] for res in restypes}
        sels = {res : uni.select_atoms(f"resname {res}", updating=True).split('residue') for res in restypes}
        for ts in uni.trajectory[frameStart:frameStop:frameStep]:
            print(ts.frame, flush=(((ts.frame-frameStart)/frameStep)%10==0))
            for res, sel in sels.items():
                #for residue in sel.split('residue'):
                for residue in sel:
                    dist_array = distances.distance_array(residue.positions, prot.positions,box=uni.dimensions)
                    addDistances[res].append(np.min(dist_array))
    
    return addDistances

def plot_hist(dists, longitudinal=False, separate=True, normalize = False, figure="distributions.png", **kwargs):
    if separate:
        fig, axis = plt.subplots(1, len(dists.keys()), figsize=(10*len(dists.keys()), 10))

        for i, res in enumerate(dists.keys()):
            axis[i].hist(dists[res], density=normalize)
            axis[i].set_title(res)
            axis[i].set_xlabel("Distance to nearest protein atom (" + r"$\mathrm{\AA}$" + ")")
            if normalize:
                axis[i].set_ylabel("Probabilty density")
            else:
                axis[i].set_ylabel("Count")
    else:
        fig, axis = plt.subplots()
        axis.set_title(" ".join(dists.keys()) + " Combined")
        axis.set_xlabel("Distance to nearest protein atom (" + r"$\mathrm{\AA}$" + ")")
        if normalize:
            axis.set_ylabel("Probabilty density")
        else:
            axis.set_ylabel("Count")

        for res in dists:
            axis.hist(dists[res], alpha=0.3, density=normalize, label = res, bins = 30)
        plt.legend()

    plt.savefig(figure)


def plot_mean(dists, figure="distributions.png", **kwargs):
    fig, axis = plt.subplots()
    items = [(k, v) for k,v in dists.items()]

    res = items[0][0].split()[0]
    axis.set_title(f"{res} Mean Distance")
    axis.set_xlabel("Frame")
    axis.set_ylabel("Mean distance to nearest protein atom (" + r"$\mathrm{\AA}$" + ")")
    
    items = list(map(lambda pair: (int(pair[0].split()[2]), np.mean(pair[1])), items))
    frame, mean = zip(*items)
    # for res in dists:
    #     axis.hist(dists[res], alpha=0.5, density=normalize, label = res, bins = 30)
    # plt.legend()
    axis.plot(frame, mean)

    plt.savefig(figure)


if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    main(ap.parse_args().infile)
