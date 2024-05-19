#!/usr/bin/env python
# -*- coding: utf-8 -*-
#"""
#Created on Tue Jan 15 11:06:39 2019
#@author: sergio Martinez
#"""
###############################################################################
# MODULES TO IMPORT -----------------------------------------------------------
###############################################################################
import os
import sys
import time
import glob
import numpy as np
import cantera as ct
from InputReader import KindOfFiles
import matplotlib.pyplot as plt
###############################################################################
# READING MECHANISMS FOLDER ---------------------------------------------------
###############################################################################
ct.suppress_thermo_warnings();
def ROP(mechanism):
    mechanism2 = mechanism.split("\\")[-1]; 
    mechanism3 = mechanism2.split(".cti")[0]
    gas = ct.Solution(mechanism)
    InputFileName = (os.getcwd()+'\\Input.inp');
    WHOLE_X,Ti,Pi,Elements,Percent,onedict = ( 
               KindOfFiles.InputFileReader(InputFileName))
    for r in range(len(Elements)):
        ELE = Elements[r]
        print("\t>Doing Reaction Path Diagram for: ",ELE)
        for m in range(len(Ti)):
           gas.TPX = Ti[m], Pi[m]*ct.one_atm, str(WHOLE_X)
           r = ct.IdealGasReactor(gas)
           sim = ct.ReactorNet([r])
           t = 0.0
           Specie = str(onedict["FUEL"])
           percentage = Percent[m]
           perper = percentage*100
           FuelConsumsion = float((r.thermo[Specie].X[0])*percentage) 
           StopFuelConsumsion =float((r.thermo[Specie].X[0]) - FuelConsumsion); 
           Fuel = r.thermo[Specie].X[0]; 
           print("\t>Fuel X initial: ",Fuel," ","100%");
           print("\t>Fuel consumption percentage criteria: ",
                 percentage*100,"%")
           print("\t>Which means, simulation will stop when Fuel X is: ",
                 StopFuelConsumsion , "reached")
           TIMES=[]; TEMP=[];PRESS=[];XFUEL=[];
           while r.thermo[Specie].X[0] > StopFuelConsumsion: 
               t += float(onedict["DELTA"])
               sim.advance(t)
               TIMES.append(sim.time)
               TEMP.append(r.T)
               PRESS.append(r.thermo.P/1E5)
               XFUEL.append(r.thermo["O2"].X[0])
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # IF YOU WANT TO PRINT ON SCREEN THE TERMODINAMIC PROPERTIES USE THE NEXT LINE
    # print(sim.time," ",gas.T," ",r.thermo['C2H4'].X[0]," ",r.thermo['CH4'].X[0])
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           element = str(ELE)
           diagram = ct.ReactionPathDiagram(gas, element)
           diagram.title ='Reaction path diagram following {0}'.format(element)
           diagram.scale=-1
           diagram.threshold = float(onedict["THRES"])
           diagram.label_threshold = float(onedict["LABEL"])
           diagram.show_details = True
           outputdata = os.getcwd()+('\\OutputData_ROP\\')
           dot_file = outputdata+str(mechanism3)+"_"+str(Ti[m])+"K_"+str(
                   onedict["FUEL"])+"_"+element+"_"+str(perper)+'%_path.dot'
           img_file = outputdata+str(mechanism3)+"_"+str(Ti[m])+"K_"+str(
                   onedict["FUEL"])+"_"+element+"_"+str(perper)+'%_path.png'
           img_path = os.path.join(os.getcwd(), img_file)
           diagram.write_dot(dot_file)
    #      print(diagram.get_data())
           print("="*90,'\n');
           print("\t>Wrote graphviz input file to '{0}'.".format(
                   os.path.join(os.getcwd(), dot_file)))
           os.system('dot {0} -Tpng -o{1} -Gdpi=200'.format(
                   dot_file, img_file))
           print("\t>Wrote graphviz output file to '{0}'.".format(
                   img_path))
           print("="*90,'\n');
           dTdt = np.append(np.diff(TEMP)/np.diff(TIMES),0)
           idt = TIMES[np.argmax(dTdt)]
           print("IDT = ",idt); print("Time = ",TIMES[-1])
           plt.plot(TIMES,XFUEL)
           #plt.show()
#+++++++++++++++++++++++++
if __name__ == "__main__":
   print("="*90,'\n');
   print("\t>RUNNING RATE OF PRODUCTION (ROP) CODE")
   print("="*90,'\n');
   mech_list = []; 
   fname = glob.glob(os.getcwd()+'\\mechanisms\\*.*');
   for i,mechanism in enumerate(fname):
       print("="*90); 
       print('\t>Running with mechanism: \n\t>',mechanism);
       print("="*90,'\n');
       ROP(mechanism)
