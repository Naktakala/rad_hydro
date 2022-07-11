import os
import subprocess

all_files = os.listdir()

time_dump_files = []
for filename in all_files:
    if filename.find("_t") > 0 and filename.find(".txt") > 0:
        time_dump_files.append(filename)

time_dump_files.sort()
for filename in time_dump_files:
    L = len(filename)
    print("Processing "+filename)
    os.system("cp " + filename + " ZRawOutput.txt")
    os.system("python3 Y0PlotRawOutput.py temperature radT")
    os.system("cp ZRawOutput.png " + filename[0:(L-4)] + ".png")