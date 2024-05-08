# -*- coding: utf-8 -*-
#"""
#Created on 2018
#@author: Sergio Martinez
# This code snippet reads voltage converted into pressure experimental
# data. It transfers it into volume data to be used as an RCM input
# volume history file to run simulations using Cantera tools.
# Important Note: please be aware that this snippet requires 
# a third software executable (segmentfit.exe). To obtain this, please
# email directly to "sergioesmartines@gmail" and place on the subject
# "software to press2volume.py" title. Thanks.
#"""
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# IMPORT MODULES                                                              +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import os
import glob
import shutil
import numpy      as np
import pandas     as pd
import cantera    as ct
import subprocess as sp
from   shutil     import copyfile
from   scipy      import optimize
import matplotlib.pyplot as plt

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CONVERTING PRESSURE TRACE TO VOLUME TRACE                                   +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def P2V(path,mech):
#    path= os.getcwd()+'\\VPRO\\'; 
    ct.suppress_thermo_warnings(); 
    names= glob.glob(path+'\\*.txt_00100_fit.txt');
    for m in range(len(names)):
        Mod=names[m].rsplit("\\", 1); 
        FileName = str(Mod[-1]);
        s= Mod[1]; sm=s.rsplit('.txt_00100_fit',1); s1=sm[0]; 
        ss = s1.split("_"); NameTest = ss[0]; Ti = ss[1];Pi=ss[2];
        PC = ss[-1] ; Pinitial = float(Pi)* 0.000986923;
        size = len(ss)
        number = int((size - 4)/2)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## MAKING NONREACTIVE AND REACTIVE LISTS AND FILLING THEM            
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        specie = []; x = []; ALL = [];REAC = []; NREAC =[];pNREAC = [];
        for n in range(number):
            specie.append(ss[int(2*n+3)])
            x.append(ss[int(2*n+4)])
        INDEX_N2 = specie.index("N2"); 
        xN2 = x[INDEX_N2]; 
        xOX = x[-1]; 
        for i in range(number):
            ALL.append(str(specie[i]) + " " + (x[i]))
        INDEX_OX = ALL.index("O2 0.0")
        NREAC = ALL.copy()
        NREAC.remove(NREAC[INDEX_OX])
        NREAC.remove(NREAC[-1])
        NEW_N2 = "N2" + " " + str(float(xN2) - float(xOX))
        pNREAC = ALL.copy()
        pNREAC.remove(NREAC[INDEX_OX])
        pNREAC.remove("O2 0.0")
        REAC = pNREAC.copy()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
######## DOING SOME MATHS FOR THE REACTIVE MOLE FRACTIONS                     
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# WRITTING THE INPUT FILE TO FEED adip2v.exe                                  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Input = NameTest+"_"+Ti+"_"+Pi+".inp"
        PressFile = NameTest+"_"+Ti+"_"+Pi+".txt"
        Output = NameTest+"_"+Ti+"_"+Pi
        copyfile(FileName,PressFile)
        output_file = open(Input, 'w')
        output_file.write("!Input file for adip2v.exe\n")
        output_file.write("!End time for the simulation\n")
        output_file.write("TIME 1\n")
        output_file.write("!Time-step for printing to diagnostic ouput file\n")
        output_file.write("DELT 1e-4\n")
        output_file.write("!Initial Temperature K\n")
        output_file.write("TEMP"+" "+Ti+"\n")
        output_file.write("!Non Reactive Mixture\n")
        for p in range(len(NREAC)):
            output_file.write("NREAC" + " " + str(NREAC[p])+"\n")
        output_file.write("!Reactive Mixture\n")
        for s in range(len(REAC)):
            output_file.write("REAC" + " " + str(REAC[s])+"\n")
        output_file.write("REAC" + " " + str(NEW_N2)+"\n")
        output_file.write("!Added to prevent convergence issues in CKPro 15101\n")
        output_file.write("REAC HE 1e-10\n")
        output_file.write("!Added to prevent convergence issues in CKPro 15101\n")
        output_file.write("REAC H2O 1e-8\n")
        output_file.close()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ProgramName = "adi_p2v.exe"
        RawFile = Input
        segment = sp.Popen([ProgramName,RawFile,str(mech),PressFile])
        wait = sp.Popen.wait(segment)
        print(wait)
        os.remove(Output+".INP_pTV.txt")
        os.remove(Output+".INP_"+str(mecha)+"_gamma.txt")
        os.remove(PressFile)
        os.remove(FileName)
        os.remove(Input)
        os.remove(Input+"_NR")
        shutil.copy(path+'\\'+Output+'.inp_r',path+'\\'+Output+'.VPRO')
#        os.rename(Input+"_R",NameTest+"_"+Ti+"_"+Pi+'.VPRO')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == "__main__":
    path_traces = os.getcwd()
    path_mechs = os.getcwd()
    mech = (glob.glob(path_mechs+"\\*.dat"));
    for i in mech:
        mecha = i.split("\\")[-1] ; 
    P2V(path_traces,mecha);
###############################################################################
