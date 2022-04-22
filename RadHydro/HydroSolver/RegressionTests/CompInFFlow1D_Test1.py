import numpy as np
import matplotlib.pyplot as plt
import typing


def ReadFile(file_name: str, offset: float = 0.0):
    file = open(file_name, "r")

    data = {}
    data["x"] = []
    data["rho"] = []
    data["p"] = []
    data["u"] = []
    data["e"] = []

    lines = file.readlines()
    for i in range(0, len(lines)):
        if i == 0:
            continue

        line = lines[i]
        words = line.split()

        if len(words) == 0:
            continue
        data["x"]  .append(float(words[1]) + offset)
        data["rho"].append(float(words[2]))
        data["p"]  .append(float(words[3]))
        data["u"]  .append(float(words[4]))
        data["e"]  .append(float(words[5]))

    file.close()
    return data

home_dir = "RadHydro/HydroSolver/RegressionTests/"

analytical = ReadFile(home_dir + "CompInFFlow1D_Test1_analytic_output.txt")
numerical  = ReadFile(home_dir + "CompInFFlow1D_Test1_output.txt", -0.5)

def PlotStuff(axs, key: str, xlim, ylim, title: str):
    axs.plot(numerical ["x"], numerical [key], label="MHM-HLLC",
                                                 linestyle="none",
                                                 marker="o",
                                                 markeredgecolor="k",
                                                 markerfacecolor="k",
                                                 fillstyle="none",
                                                 markersize=3)
    axs.plot(analytical["x"], analytical[key], label="Analytical", color="k")
    axs.legend()
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)
    axs.set_xlabel("x")
    axs.set_ylabel(title)

fig, axs = plt.subplots(2,2,figsize=[8.0,8.0],dpi=200.0)
PlotStuff(axs[0,0],"rho",[-0.5,0.5],[0.0,1.2], "Density")
PlotStuff(axs[0,1],"u",[-0.5,0.5],[0.0,1.2], "Velocity")
PlotStuff(axs[1,0],"p",[-0.5,0.5],[0.0,1.2], "Pressure")
PlotStuff(axs[1,1],"e",[-0.5,0.5],[1.0,3.4], "Energy")

plt.savefig("CompInFFlow1D_Test1_output.png")