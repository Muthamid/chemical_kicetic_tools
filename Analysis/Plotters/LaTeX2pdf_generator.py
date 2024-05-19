#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:23:58 2018
@author: Sergio Martinez
"""
import subprocess
import time
import os
import sys
import glob
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
from matplotlib import colors as mcolors
from cycler import cycler
import time
timestr = time.strftime("%d_%m_%Y")
print(timestr)
# ---------------------------------------------------------------
#from STplotter.py import 
# ---------------------------------------------------------------
CWD       = os.getcwd();
MONOS        = ["CH4","C2H4","C2H6","C3H8"]
BINARIES     = ["CH4_C2H4","CH4_C2H6","C2H4_C2H6","C2H4_C3H8","C2H6_C3H8"]
TERNARIES    = ["CH4_C2H4_C2H6"]
QUATERNARIES = ["CH4_C2H4_C2H6_C3H8"]
STreader  = pd.read_csv(CWD+"\\ST_data.LateX.out");
RCMreader = pd.read_csv(CWD+"\\RCM_data.LateX.out");
FileNameInit = "LatexOutputFile_"+timestr+".tex"
f = open(FileNameInit, "w")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\documentclass{article}                          \n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("% REQUESTED PACKAGES                              \n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\usepackage{fullpage} % USER REQUESTED PACKAGE   \n")
f.write("\\usepackage{amssymb} % USER REQUESTED PACKAGE    \n")
f.write("\\usepackage{amsmath} % USER REQUESTED PACKAGE    \n")
f.write("\\usepackage{rotating} % USER REQUESTED PACKAGE   \n")
f.write("\\usepackage{ctable} % USER REQUESTED PACKAGE     \n")
f.write("\\usepackage{longtable} % USER REQUESTED PACKAGE  \n")
f.write("\\usepackage{lscape} % USER REQUESTED PACKAGE     \n")
f.write("\\usepackage{graphicx} % USER REQUESTED PACKAGE   \n")
f.write("\\usepackage{color} % USER REQUESTED PACKAGE      \n")
f.write("\\usepackage{subfig} % USER REQUESTED PACKAGE     \n")
f.write("\\usepackage{caption} % USER REQUESTED PACKAGE    \n")
f.write("\\usepackage{hyperref} % USER REQUESTED PACKAGE   \n")
f.write("\\usepackage{makeidx} % USER REQUESTED PACKAGE   \n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\makeindex\n")
f.write("%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("% TITLE INFORMATION                                           \n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\title{\\textbf{\\break \\Huge Combustion Chemistry Centre}    \n")   
f.write("\\break \\textbf{\LARGE National University of Ireland Galway}\n")
f.write("\\\\")
f.write("\\textbf{NUIGMech1.2\_C4 \\\ VS \\\ NUIGMech1.3\_C4 \\\ Mechanism \\\ validation}\n")
f.write("\\break \\textit{Chemical Kinetics Report:\n")
f.write("\\\ Single C1, C2, and C3\n")
f.write("\\\ Binary, Ternary and Quaternary blends\n")
f.write("\\\ of gaseous hydrocarbons}}\n")
f.write("\\date{\\today}\n")
f.write("\\author{Sergio Martinez}\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\begin{document}\n")
f.write("\\maketitle\n")
f.write("\\tableofcontents{}\n")
#f.write("\\listoffigures\n")
f.write("\\clearpage\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("% Introduction\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
#f.write("\\section{Introduction}\n")
#f.write("\\subsection{How to cite NUIGMech1.1:}\n")
#f.write("\\paragraph{NUIGMech1.1 - October 2nd 2020 \\\ Important notice: NUIGMech1.0 Published first in citation no 1 below. Please cite all papers below when referencing it.}\n")
#f.write("\\paragraph{1. M. Baigmohammadi, V. Patel, S. Nagaraja, A. Ramalingam, S. Martinez, S. Panigrahy, A. Mohamed, K.P. Somers, U. Burke, K.A. Heufer, A. Pekalski, H.J. Curran, Comprehensive experimental and simulation study of the ignition delay time characteristics of binary blended methane, ethane, and ethylene over a wide range of temperature, pressure, equivalence ratio, and dilution, Energy Fuels 34 (2020) 8808-8823.\n")
#f.write("\\\ 2. S.S. Nagaraja, J. Liang, S. Dong, S. Panigrahy, A.B. Sahu, G. Kukkadapu, W.J. Pitz, H.J. Curran, A hierarchical single-pulse shock tube pyrolysis study of C2-C6 1-alkenes, Combustion and Flame 219 (2020) 456-466.                                                                                                                                                                                                    \n")
#f.write("\\\ 3. N. Lokachari, S. Panigrahy, G. Kukkadapu, G. Kim, T. MacDougall, S. Vasu, W.J. Pitz, H.J. Curran, The influence of isobutene kinetics on the reactivity of di-isobutylene and iso-octane, Combustion and Flame 222 (2020) 186-195.                                                                                                                                                                                      \n")
#f.write("\\\ 4. S. Panigrahy, J. Liang, S.S. Nagaraja, Z. Zuo, G. Kim, T. MacDougall, S.S. Vasu, H.J. Curran, A comprehensive experimental and improved kinetic modeling study on the pyrolysis and oxidation of propyne, Proceedings of the Combustion Institute 38 (2021) accepted.                                                                                                                                                   \n")
#f.write("\\\ 5. A. Abd El-Sabor Mohamed, S. Panigrahy, A.B. Sahu, G. Bourque, H.J. Curran, An experimental and modeling study of the auto-ignition of natural gas blends containing C1-C7 n-alkanes, Proceedings of the Combustion Institute 38 (2021) accepted.                                                                                                                                                                        \n")
#f.write("\\\ 6. S.S. Nagaraja, J. Power, G. Kukkadapu, S. Dong, S.W. Wagnon, W.J. Pitz, H.J. Curran, A single pulse shock tube study of pentene isomer pyrolysis, Proceedings of the Combustion Institute 38 (2021) accepted.                                                                                                                                                                                                           \n")
#f.write("\\\ 7. S. Dong, K. Zhang, P.K. Senecal, G. Kukkadapu, S.W. Wagnon, S. Barrett, N. Lokachari, S. Panigrahy, W.J. Pitz, H.J. Curran, A comparative reactivity study of 1-alkene fuels from ethylene to 1-heptene, Proceedings of the Combustion Institute 38 (2021) accepted.                                                                                                                                                    \n")
#f.write("\\\ 8. S. Dong, K. Zhang, E.M. Ninnemann, A. Najjar, G. Kukkadapu, J. Baker, F. Arafin, Z. Wang, W.J. Pitz, S.S. Vasu, S.M. Sarathy, P.K. Senecal, H.J. Curran, A comprehensive experimental and kinetic modeling study of 1- and 2-pentene, Combustion and Flame (2021) accepted.                                                                                                                                             }\n")
#f.write("\\subsection{Main changes: from NUIGMech1.0 to NUIGMech1.1}\n")
#f.write("\\paragraph{This is just an illustrative document to compare mechanism's performance.}\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("% Single Fuels Section:\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\section{Single Fuel Section}\n")
f.write("\\subsection{Shock Tube Validation:}\n")
f.write("\\begin{center}\n")
for k in range(len(STreader)):
       FNam  = ((str(STreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == MONOS[0]:
             f.write(str(STreader.iloc[k,0])+"\n")
       if    str(FNam) == MONOS[1]:
             f.write(str(STreader.iloc[k,0])+"\n")
       if    str(FNam) == MONOS[2]:
             f.write(str(STreader.iloc[k,0])+"\n")
       if    str(FNam) == MONOS[3]:
             f.write(str(STreader.iloc[k,0])+"\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\subsection{Rapid Compression Machine Validation:}\n")
for k in range(len(RCMreader)):
       FNam  = ((str(RCMreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == MONOS[0]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
       if    str(FNam) == MONOS[1]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
       if    str(FNam) == MONOS[2]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
       if    str(FNam) == MONOS[3]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
f.write("\\end{center}\n")
f.write("\\clearpage\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("% Binary Blends Section:\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\section{Binary Blends Section}\n")
f.write("\\subsection{Shock Tube Validation:}\n")
f.write("\\begin{center}\n")
for k in range(len(STreader)):
       FNam  = ((str(STreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == BINARIES[0]:
             f.write(str(STreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[1]:
             f.write(str(STreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[2]:
             f.write(str(STreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[3]:
             f.write(str(STreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[4]:
             f.write(str(STreader.iloc[k,0])+"\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\subsection{Rapid Compression Machine Validation:}\n")
for k in range(len(RCMreader)):
       FNam  = ((str(RCMreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == BINARIES[0]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[1]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[2]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[3]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
       if    str(FNam) == BINARIES[4]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
f.write("\\end{center}\n")
f.write("\\clearpage\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("% Ternary Blends Section:\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\section{Ternary Blends Section}\n")
f.write("\\subsection{Shock Tube Validation:}\n")
f.write("\\begin{center}\n")
for k in range(len(STreader)):
       FNam  = ((str(STreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == TERNARIES[0]:
             f.write(str(STreader.iloc[k,0])+"\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\subsection{Rapid Compression Machine Validation:}\n")
for k in range(len(RCMreader)):
       FNam  = ((str(RCMreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == TERNARIES[0]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
f.write("\\end{center}\n")
f.write("\\clearpage\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("% Quaternary Blends Section:\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\section{Quaternary Blends Section}\n")
f.write("\\subsection{Shock Tube Validation:}\n")
f.write("\\begin{center}\n")
for k in range(len(STreader)):
       FNam  = ((str(STreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == QUATERNARIES[0]:
             f.write(str(STreader.iloc[k,0])+"\n")
f.write("%+++++++++++++++++++++++++++++++++++++++++++++++++\n")
f.write("\\subsection{Rapid Compression Machine Validation:}\n")
for k in range(len(RCMreader)):
       FNam  = ((str(RCMreader.iloc[k,0]).split('/'))[-4]);
       if    str(FNam) == QUATERNARIES[0]:
             f.write(str(RCMreader.iloc[k,0])+"\n")
f.write("\\end{center}\n")
f.write("\\clearpage\n")
f.write("\\printindex{}\n")
f.write("\\end{document}")
f.close()
for k in range(2):
    subprocess.call(["pdflatex",FileNameInit])
#
#time.sleep(5)
#
#pdfreader = "C:\\Program Files (x86)\\Adobe\\Acrobat Reader DC\\Reader\\AcroRd32.exe"
#pre, ext = os.path.splitext(FileNameInit) 
#pdfoutput = pre + ".pdf"
#
#subprocess.call([pdfreader,pdfoutput])
