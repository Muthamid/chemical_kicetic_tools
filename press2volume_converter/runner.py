# -*- coding: utf-8 -*-
#"""
#Created in 2018
#@author: Sergio Martinez
# This code snippet reads voltage converted into pressure experimental
# data. It transfers it into volume data to be used as an RCM input
# volume history file to run simulations using Cantera tools.
# Important Note: please be aware that this snippet requires 
# a third software executable (segmentfit.exe). To obtain this, please
# email directly to "sergioesmartines@gmail" and place on the subject
# "software to press2volume.py" title. Thanks.
#"""
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# IMPORT MODULES
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import os
from scipy import optimize
import numpy as np
import pandas as pd
import cantera as ct
import subprocess as sp
import glob
import matplotlib.pyplot as plt
np.set_printoptions(precision=5)
ct.suppress_thermo_warnings(); ReadPath = os.getcwd()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CONVERTING PRESSURE TRACE TO VOLUME TRACE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def P2V(path):
#    path= os.getcwd()+'\\VPRO\\'; 
    names= glob.glob(path+'\\*.txt');
    for m in range(len(names)):
        Mod=names[m].rsplit("\\", 1)
        FileName = str(Mod[-1])
        s= Mod[1]; sm=s.rsplit('.csv.txt',1); s1=sm[0]; 
        ss = s1.split("_"); Pi=ss[2]; Pinitial = float(Pi)/1000;
        ProgramName = "segmentfit.exe"
        FileName = names[m]; NUMBER = "00100";
        segment = sp.Popen([ProgramName, FileName,str(NUMBER)])
        wait = sp.Popen.wait(segment)
        print(wait)
        FitFile = FileName+"_"+str(NUMBER)+"_fit.txt"
        os.remove(FileName+"_00100_fit0x.txt")
        os.remove(FileName+"_00100_fit0.txt")
        os.remove(FileName+"_00100_fitd.txt")
        os.remove(FileName+"_00100_fitx.txt")
        os.remove(FileName+"_error.txt")
        file = pd.read_table(FitFile,names=['col'], skipinitialspace=True)
        df = pd.DataFrame(file.col.str.split('  ',1).tolist(),columns = ['time','press'])
        df["press"] = df.press.apply(lambda x: ( float(x)))
        D = Pinitial - df["press"][0]
        df["press"] = df.press.apply(lambda x: ((float(x) + D)*0.986923))
        df.to_csv(FitFile,index=False,header=False,sep="\t")
        print(df["press"])
if __name__ == "__main__":
    path = os.getcwd()
    P2V(path);
########################################################################
