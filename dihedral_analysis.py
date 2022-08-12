import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
from myTools.misc import readSimpleInput
from pickAQuote import pickAQuote
from MDAnalysisTests.datafiles import GRO, XTC
from matplotlib.patches import Rectangle
from time import perf_counter_ns


def dihedral_analysis(infile):
    print("\nProcessing input...")
    args = readSimpleInput(infile)

    required = ["top", "traj", "frameStep", "out", "plot"]

    # check for a bunch of possible errors
    errors = []
    for item in required:
        try:
            vars(args)[item]
        except AttributeError:
            errors.append(f"ERROR: missing argument '{item}'")

     # kill warning about unknown atom masses, useless
    mda.warnings.simplefilter("ignore")
    if args.top.endswith(".top"):
        try:
            print("reading GROMACS topology...")
            uni = mda.Universe(args.top, args.traj, topology_format="ITP")
        except Exception:
            print("GROMACS topology: make sure you have #included full paths to force field in topology!")
            raise RuntimeError(".top files are read as GROMACS topologies, change topology extension")
    else:
        print("reading structure/topology...")
        uni = mda.Universe(args.top, args.traj)

    # uni = mda.Universe(GRO, XTC)

    prot = uni.select_atoms("protein")

    frames = len(uni.trajectory)//args.frameStep

    fig, ax = plt.subplots()
    print("Begin Ramachandran Analysis", flush=True)
    start = perf_counter_ns()
    r = Ramachandran(prot)
    output = r.run(step = args.frameStep)
    ax = r.plot(ax, s = 0.1)
    end = perf_counter_ns()
    print(f"{frames} frames analyzed in {(end - start)/10**9} s", flush=True)


    # Beta Sheet Criteria
    # https://journals.iucr.org/d/issues/2002/05/00/gr2189/index.html
    # -180 < phi < -45, 45 < psi < 225

    # count = len(output.results.angles[0])
    # print_str = ""
    # with open(args.out, "w+") as o:
    #     start_write = perf_counter_ns()
    #     for frame in output.results.angles:
    #         angle_list = []
    #         beta = 0
            
    #         for phi, psi in frame:
    #             angle_list.append(f"{phi},{psi}")
    #             if -180 < phi < -45 and 45 < psi < 225:
    #                 beta += 1
    #         print_str+=f"{beta/count*100:.1f}% Beta Sheet\n"
    #         o_str = " ".join(angle_list)
    #         o.write(o_str+"\n")
    # print(print_str)
    # print(f"{frames} frames written in {(perf_counter_ns() - start_write)/10**9} s")

    count = len(output.results.angles[0])*len(output.results.angles)
    beta = 0
    for frame in output.results.angles:
        for phi, psi in frame:
            if -180 < phi < -45 and (psi > 45 or psi < -135):
                beta += 1
    print(f"{beta/count*100:.1f}% Beta Sheet\n")

    ax.set_title(" ".join(args.plot.split(".")[0].split("_")))
    ax.add_patch(Rectangle((-180,45),135,135, edgecolor = "red", fill = False))
    ax.add_patch(Rectangle((-180,-180),135,45, edgecolor = "red", fill = False))
    ax.text(60, 150, f"{beta/count*100:.1f}% Beta Sheet", fontdict={"color": "red"})
    plt.savefig(args.plot)

    pickAQuote("./quotes.txt")

    

if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    dihedral_analysis(ap.parse_args().infile)