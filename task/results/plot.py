import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import re;
import os, sys

CONTENT_REGEX = re.compile("Dimentions: (?P<size>\d+) x.*\niteration (?P<iter>\d+) .*\n\s+(?P<time>\d+\.\d+)s", re.M)
FILENAME_REGEX = re.compile("(?P<node>\d+)\.ppch_(?P<type>(?:[a-z]+_)?[a-z]+).*?\.\d+\.stdout")

ALLOWED_DIRS  = ["by_size", "by_node"]
ALLOWED_TYPES = ["mpi", "mpi_omp"]

def get_dir_name(path):
    return path[::-1].split(os.sep, 1)[0][::-1]

def parse_results_tree():
    res = {}
    for t in ALLOWED_TYPES:
        res[t] = {}
        for d in ALLOWED_DIRS:
            res[t][d] = []

    for root, _, files in os.walk(os.getcwd()):
        dirname = get_dir_name(root)
        if dirname not in ALLOWED_DIRS:
            continue
        #print("{0} has dirname {1}".format(root, dirname))

        for f in files:
            fmatch = FILENAME_REGEX.match(f)
            if not fmatch:
                continue
            #print("{0} matches regex!".format(f))

            fd = fmatch.groupdict()
            nd, ft = fd["node"], fd["type"]
            #print("nodeid {0}, type {1}".format(nd, ft))

            fh = open(os.path.join(root, f), "r")
            data = fh.read()
            mdata = CONTENT_REGEX.match(data)
            if not mdata:
                print("WRONG FILE SYNTAX: {0}/{1}", root, f)
            mdatad = mdata.groupdict()

            res[ft][dirname] += [{"node": int(nd), "size": int(mdatad["size"]), "iter": int(mdatad["iter"]), "time": float(mdatad["time"])}]
    return res

def plot():
    matplotlib.rcParams.update({'font.size': 22})
    fig, plts = plt.subplots(2, 1)

    plotdata = parse_results_tree()
    #print(plotdata)

    for i, testt in enumerate(ALLOWED_DIRS):
        plts[i].set_title(testt)
        plts[i].set_ylabel("Time (s)")
        
        mpi_legend = mpatches.Patch(color='red', label='MPI')
        mpi_omp_legend = mpatches.Patch(color='green', label='MPI + OMP')
        plts[i].legend(handles=[mpi_legend, mpi_omp_legend])

        for mpit in ALLOWED_TYPES:
            color = "r" if mpit == "mpi" else "g"
            if testt == "by_size":
                entries = sorted(plotdata[mpit][testt], key = lambda x: x["size"])
                Y = [e["time"] for e in entries]
                X = [e["size"] for e in entries]
                plts[i].set_xlabel("Matrix size")
                plts[i].plot(X, Y, color = color, linewidth = 2)
            else:
                entries = sorted(plotdata[mpit][testt], key = lambda x: x["node"])
                Y = [e["time"] for e in entries]
                X = [e["node"] for e in entries]
                print("X = ", X, "Y = ", Y)
                plts[i].set_xlabel("MPI node count")
                plts[i].plot(X, Y, color = color, linewidth = 2)
    plt.show()

def iter_plot():
    matplotlib.rcParams.update({'font.size': 22})
    plotdata = parse_results_tree()

    entries = sorted(plotdata["mpi"]["by_size"], key = lambda x: x["size"])
    Y = np.array([e["iter"] for e in entries])
    X = np.array([e["size"] for e in entries])
    plt.xlabel("Matrix size")
    plt.ylabel("Iterations")
    plt.title("Iterations by size")
    plt.plot(X, Y, "b", linewidth = 2)

    #k = 0.25
    #SX = np.power(X - X[0], 2) * k + Y[0]
    #plt.plot(X, SX, "r--", linewidth = 2)

    iter_legend = mpatches.Patch(color='blue', label='iter f()')
    sugg_legend = mpatches.Patch(color='red', label='linear(size^2)')
    plt.legend(handles=[iter_legend])#, sugg_legend])

    plt.show()