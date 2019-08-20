#!/usr/bin/env python3
import re
import os, os.path
import shutil
import glob
import subprocess
import multiprocessing
from multiprocessing import Pool
import warnings



# color definition
class TColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def convert_pdf(infile, outdir="./img"):
    # Convert svg to pdf
    # print(infile)
    base_name = os.path.basename(infile)
    #FIRST round svg --> ps
    if ("tmp" in base_name) or (".svg" not in base_name):      # temporary svg files
        print(TColors.OKBLUE +
              "File {} ommited.".format(infile)
              +TColors.ENDC)
        return
    outfile = os.path.join(os.path.abspath(outdir), base_name)
    # outfile = outfile.replace(".svg", ".pdf")
    outfile = outfile.replace(".svg", ".png")
    infile = os.path.abspath(infile)
    #outfile = outfile.replace(".svg", ".png")
    print(outfile)
    program = "/Applications/Inkscape.app/Contents/Resources/bin/inkscape"
    params = ["--without-gui", "--export-area-page",
    ]
    io = ["--file={}".format(infile),
          "--export-dpi=600",
          "--export-png={}".format(outfile)]
    success = subprocess.call([program, *params, *io])
    if success != 0:
        warnings.warn(TColors.FAIL + "File {} cannot be converted!".format(infile) + TColors.ENDC)
    else:
        print(TColors.OKGREEN +
              "File {} converted to ps on thread {}.".format(infile, multiprocessing.current_process())
              +TColors.ENDC)

    #SECOND round ps to pdf using ps2pdf
    program = "convert"
    infile = outfile
    outfile = infile.replace(".png", ".pdf")
    io = [infile, outfile]
    success = subprocess.call([program, *io])
    if success != 0:
        warnings.warn(TColors.FAIL + "File {} cannot be converted!".format(infile) + TColors.ENDC)
    else:
        print(TColors.OKGREEN +
              "File {} converted to pdf on thread {}.".format(infile, multiprocessing.current_process())
              +TColors.ENDC)

# Convert all the pdf files
cur_path = os.path.dirname(__file__)
img_path = os.path.join(cur_path, "./img")

if __name__ == "__main__":
    file_list = []
    # TeX process
    # for ifile in glob.glob("*.tex"):
        # tex_process(ifile)
        # print(TColors.OKBLUE + "Converted TeX file: {}".format(ifile) + TColors.ENDC)

    # PDF Process 
    for ifile in glob.glob(os.path.join(img_path, "*.svg")):
        this_img_dir = os.path.dirname(ifile)
        file_list.append((ifile, this_img_dir))
    N_cores = multiprocessing.cpu_count()
    # multicore
    with Pool(N_cores) as p:
        p.starmap(convert_pdf, file_list)
