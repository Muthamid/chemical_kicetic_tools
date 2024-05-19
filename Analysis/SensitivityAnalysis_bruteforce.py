#!/usr/bin/env python
# -*- coding: utf-8 -*-
#"""
#Created on Tue Jan 15 11:06:39 2019
#@author: sergio Martinez
#"""
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# IMPORTING MODULES - LET'S ROCK!
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import os
import numpy as np
import cantera as ct
import pandas as pd
import glob
from InputReader import KindOfFiles
import itertools
import multiprocessing
from Definitions.Definitions import DefinitionStandar
import csv
from time import time
import datetime
import platform

gases = {}

def idtST(mech, T, P, X, factor):
    ct.suppress_thermo_warnings();
    mech, T, P, X, factor 
    gas = ct.Solution(mech)
    gas.TPX = T, P, X
    gas.set_multiplier(1.0)
    for i in range(gas.n_reactions):
       gas.set_multiplier(factor[i],i)
    r = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([r])
    time = []
    temp = []
    pres = []
    states = ct.SolutionArray(gas, extra=['t'])
    print("\t>Ti = {0} K, Pi= {1} atm".format(T,round(P*9.86923e-6,3)))
    while sim.time < 1.0 and r.T < (400+T): 
          sim.step()
          time.append(sim.time)
          temp.append(r.T)
          pres.append(r.thermo.P)
          states.append(r.thermo.state, t=sim.time)
    time = np.array(time)
    temp = np.array(temp)
    pres = np.array(pres)
    diff_temp = np.diff(temp)/np.diff(time)
    dpdt      = np.diff(pres)/np.diff(time)
    ign_temp   = np.argmax( diff_temp )
    ign_pres  = np.argmax( dpdt )
    ign = time[ign_pres]
    if ign == 0:
        ign = time[ign_temp]
    return (ign, T, P, X)

def idtST2(mech, T, P, X, factor):
    ct.suppress_thermo_warnings();
    mech, T, P, X, factor 
    gas = ct.Solution(mech)
    gas.TPX = T, P, X
    gas.set_multiplier(1.0)
    for i in range(gas.n_reactions):
       gas.set_multiplier(factor[i],i)
    r = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([r])
    time = []
    temp = []
    pres = []
    states = ct.SolutionArray(gas, extra=['t'])
    #print("\t>Ti = {0} K, Pi= {1} atm".format(T,P))
    while sim.time < 1.0 and r.T < (400+T): 
          sim.step()
          time.append(sim.time)
          temp.append(r.T)
          pres.append(r.thermo.P)
          states.append(r.thermo.state, t=sim.time)
    time = np.array(time)
    temp = np.array(temp)
    pres = np.array(pres)
    diff_temp = np.diff(temp)/np.diff(time)
    dpdt      = np.diff(pres)/np.diff(time)
    ign_temp   = np.argmax( diff_temp )
    ign_pres  = np.argmax( dpdt )
    ign = time[ign_pres]
    if ign == 0:
        ign = time[ign_temp]
    return (ign)

def ign_RCM(mech, factor):	
      ct.suppress_thermo_warnings()
      gato = '#' * 90
      DefinitionStandar.SOME_STUFF(90, 'cti', 'xml')
      path = os.getcwd() + '\\FILES\\'
      names = glob.glob(path + '*.inp_r')
      Ti, Pi, fuelname, whole_X, volname, TMP = (DefinitionStandar.
                                                   FileFormat(names, path))
      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #for i in range(len(names)):
      file_name = TMP[0]
      df = pd.read_csv(file_name)
      vpro_time = pd.to_numeric(df.iloc[:, 0])
      gas = ct.Solution(mech)
      gas.set_multiplier(1.0)
      for i in range(gas.n_reactions):
          gas.set_multiplier(factor[i],i)
      gas.TPX = Ti[0], Pi[0], whole_X[0]
      print('Fuel and mole fractions:\t\n', whole_X[0])
      print(gato)
      print(gas())
      print(gato)
      print("""RUNNING NOW:\nTi[K]={:<10.2f}\nPi[mbar]={:<10.2f}
              \nVprofile NAME={:<10s}\n..."""
                  .format(Ti[0], Pi[0] / 1E2, volname[0]))
      print(gato)
      r = ct.IdealGasReactor(gas)
      env = ct.Reservoir(ct.Solution('air.xml'))
      w = ct.Wall(r, env)
      w.set_velocity(DefinitionStandar.variable_volume_velocity(df))
      sim = ct.ReactorNet([r])
      end_time = float(vpro_time.iloc[-1] * 10)
      Deltat   = float(vpro_time.iloc[1] / 100)
      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      TIMES = []
      TEMPERATURES = []
      PRESSURES = []
      V = [] 
      t = 0
      while t < (end_time):
            t += Deltat
            sim.advance(t)
            TIMES.append(sim.time)
            TEMPERATURES.append(gas.T)
            PRESSURES.append(gas.P / 1E5)
            V.append(r.volume)
      EocArg = np.argmin(V)
      EocTime = TIMES[EocArg]
      dpdt = np.append(np.diff(PRESSURES) / np.diff(TIMES), 0)
      dTdt = np.append(np.diff(TEMPERATURES) / np.diff(TIMES), 0)
      eoc = EocTime
      IdtDpdt = TIMES[np.argmax(dpdt)] - eoc
      if IdtDpdt < 0:
            IdtDpdt = 0.0
      IdtDTdt = TIMES[np.argmax(dTdt)] - eoc
      if IdtDTdt < 0:
            IdtDTdt = 0.0
      gas.set_multiplier(1.0) 
      return (IdtDpdt,Ti[0], Pi[0], whole_X[0])

def ign_RCM2(mech, T, P, X, factor):	
      ct.suppress_thermo_warnings()
      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #for i in range(len(names)):
      CWD = os.getcwd()
      TMP = glob.glob(CWD+"\\FILES\\*.tmp")
      file_name = TMP[0]
      df = pd.read_csv(file_name)
      vpro_time = pd.to_numeric(df.iloc[:, 0])
      gas = ct.Solution(mech)
      gas.set_multiplier(1.0)
      for i in range(gas.n_reactions):
          gas.set_multiplier(factor[i],i)
      gas.TPX = T, P, X
      print('Fuel and mole fractions:\t\n', X)
      print(gas())
      #print("""RUNNING NOW:\nTi[K]={:<10.2f}\nPi[mbar]={:<10.2f}\n..."""
      #            .format(T, P/1E2))

      r = ct.IdealGasReactor(gas)
      env = ct.Reservoir(ct.Solution('air.xml'))
      w = ct.Wall(r, env)
      w.set_velocity(DefinitionStandar.variable_volume_velocity(df))
      sim = ct.ReactorNet([r])
      end_time = float(vpro_time.iloc[-1] * 10)
      Deltat = float(vpro_time.iloc[1] / 100)
      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      TIMES = []
      TEMPERATURES = []
      PRESSURES = []
      V = [] 
      t = 0
      while t < (end_time):
            t += Deltat
            sim.advance(t)
            TIMES.append(sim.time)
            TEMPERATURES.append(gas.T)
            PRESSURES.append(gas.P / 1E5)
            V.append(r.volume)
      EocArg = np.argmin(V)
      EocTime = TIMES[EocArg]
      dpdt = np.append(np.diff(PRESSURES) / np.diff(TIMES), 0)
      dTdt = np.append(np.diff(TEMPERATURES) / np.diff(TIMES), 0)
      eoc = EocTime
      IdtDpdt = TIMES[np.argmax(dpdt)] - eoc
      if IdtDpdt < 0:
            IdtDpdt = 0.0
      IdtDTdt = TIMES[np.argmax(dTdt)] - eoc
      if IdtDTdt < 0:
            IdtDTdt = 0.0
      gas.set_multiplier(1.0) 
      return (IdtDpdt)

def BF(args):
    mech, k, dk, factor, T, P, X, ign0, FileInputName, FacilityName = args
    gas = ct.Solution(mech)
    gas.TPX = T, P, X; 
    reactions = [];   m = gas.n_reactions; 
    ds = pd.DataFrame(data=[], index=gas.reaction_equations(range(m)))
    ds["index"] = ""
    ds["bruteforce"] = ""
    all_reactions = ct.Reaction.listFromFile(mech)
    for R in all_reactions:
        reactions.append(R)
    factor[k]  = 1+dk
    if   "ST" in FacilityName:
         ign = idtST2(mech,T,P,X,factor)
    elif "RCM" in FacilityName:
         ign = ign_RCM2(mech, T, P, X, factor)
    factor[k]  = 1.0
    ds["bruteforce"][k] = (ign-ign0)/(ign0*dk) 
    ds["index"][k]      = k
    SIDTs=open(str(FileInputName)+".txt","a")
    SIDTs.writelines(str(reactions[k]))
    SIDTs.writelines('	')
    SIDTs.writelines(str(ds["index"][k]))
    SIDTs.writelines('	')
    SIDTs.writelines(str(ds["bruteforce"][k])+'\n')
    SIDTs.close()

def init_process(mech):
    """
    This function is called once for each process in the Pool. We use it to
    initialize any Cantera objects we need to use.
    """
    ct.suppress_thermo_warnings();
    gases[mech] = ct.Solution(mech)

def parallel(mech, predicate, nProcs, nTemps,dk,factor,T,P,X,ign0,FileInputName,FacilityName): 
    """
    Call the function ``predicate`` on ``nProcs`` processors for ``nTemps``
    different temperatures.
    """
    ct.suppress_thermo_warnings();
    pool = multiprocessing.Pool(processes=nProcs,
                                initializer=init_process,
                                initargs=(mech,))
    reactionslist = range(nTemps);
    y = pool.map(predicate,
                 zip(itertools.repeat(mech),
                     reactionslist,
                     itertools.repeat(dk),
                     itertools.repeat(factor),
                     itertools.repeat(T),
                     itertools.repeat(P),
                     itertools.repeat(X),
                     itertools.repeat(ign0),
                     itertools.repeat(FileInputName),
                     itertools.repeat(FacilityName)))
    #pool.close()
    #pool.join()
    #pool.terminate()
    return (y)

def SENS(mechanism):
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# RCM DEFINITION: IDT CALCULATION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++
# RUNNING THE WHOLE SENSITIVITY ANALYSIS
#+++++++++++++++++++++++++++++++++++++++
# Global storage for Cantera Solution objects
 def SecondRound(m,gas, factor, dk, reactions, nProcs,T,P,X,mechanism,ds,FileInputName,FacilityName):
    if   "ST" in FacilityName:
         ign0, TT, PP, XX = idtST(mechanism,T,P,X,factor)
    elif "RCM" in FacilityName:
         ign0, TT, PP, XX = ign_RCM(mechanism, factor)
    print("Ignition Delay is: {:.4f} ms".format(ign0*1000))
    factor = np.ones( ( gas.n_reactions, 1 ) )
    print('Start Brute Force')
    print("\t> Running with {0} proccessors ".format(nProcs))
    ds = parallel(mechanism, BF, nProcs, m, dk,factor,TT,PP,XX,ign0,FileInputName,FacilityName)
    return (ds)

 # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 # CODE TO RUN: ALL SYSTEM'S FEATURES
 # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ct.suppress_thermo_warnings();
 gas = ct.Solution(mechanism)
 all_reactions = ct.Reaction.listFromFile(mechanism)
 reactions = []
 for R in all_reactions:
    reactions.append(R)
 InputFileName = (os.getcwd()+'\\Input.inp');
 (WHOLE_XF,WHOLE_XO,Ti,Pi,DK,onedict) = (KindOfFiles.InputFileReader(InputFileName))
 for mm in range(len(Ti)):
  temp = Ti[mm]
  pres = Pi[mm]
  fuel = str(WHOLE_XF)
  oxidizer = str(WHOLE_XO)
  simtype = 'SV'
  mecha1  = str(mechanism).split("\\")[-1]
  mecha2  = str(mecha1).split(".cti")[0]
  Path   = os.getcwd()
  Output = Path +"\\"+ str(mecha2)
  if not os.path.isdir(Output):
         os.makedirs(Output)
  FacilityName  = str(onedict["FACILITY"])+"_Ti_"+str(temp)+"_pi_"+str(pres)
  FileInputName = Output + "\\" + FacilityName+"_"+str(onedict["NAMEFILE"])+"_"+ str(mecha2)
  dk = DK[mm]
  X  = fuel+','+oxidizer
  P  = pres*ct.one_atm
  T  = temp
  gas.TPX = T, P, X
  gas.equilibrate(simtype)
  ds = pd.DataFrame(data=[], index=gas.reaction_equations(range(10)))
  m = gas.n_reactions
  cpus = int(multiprocessing.cpu_count())
  if   "auto" in onedict["WORKERS"]:
       nProcs = cpus
  elif "semi" in onedict["WORKERS"]:
       nProcs = cpus-1
  else:
       nProcs = int(onedict["WORKERS"])
  pd.options.display.float_format = '{:,.2e}'.format
  ds["index"] = ""
  ds["bruteforce"] = ""
  factor = np.ones( gas.n_reactions )
  SecondRound(m,gas,factor,dk,reactions,nProcs,T,P,X,mechanism,ds,FileInputName,FacilityName)
  converter = pd.read_csv(str(FileInputName)+".txt",sep="\t",names=["reactions","index","bruteforce","bruteforcetabs"])
  df        = pd.DataFrame(data=converter)
  df.to_csv(str(FileInputName)+".csv",index=False)
 return (str(FileInputName),m,nProcs)
#++++++++++++++++++++++++++++++
if __name__ == "__main__":
   print("="*90,'\n');
   print("\t> RUNNING SENSITIVITY ANALYSIS MULTIPROCESSING\n\t> Update:  30/05/2019. S.M.& C3-team NUIGalway,IE")
   print("="*90,'\n');

   mech_list = []; 
   fname = glob.glob(os.getcwd()+'\\Mechanisms\\*.*');
   for i,mechanism in enumerate(fname):
        t1    = time()
        print("="*90); 
        print('\t> Running with mechanism: \n\t>',mechanism);
        print("="*90,'\n');
        outname,numreacs,nProcs = SENS(mechanism)
        t2       = time()
        clock    = t2-t1
        today    = datetime.date.today()
        mecha1   = str(mechanism).split("\\")[-1]
        mecha2   = str(mecha1).split(".cti")[0]
        machina  = platform.machine()
        version  = platform.version()
        uname    = platform.uname()
        system   = platform.system()
        procs    = nProcs
        platformA = platform.platform()
        if   clock <= 60:
             reloj = round(clock,2)
             print("\t> Simulation took:{0:0.2f} seconds".format(reloj))
             print("Printing log.txt file...")
             SIDTs=open(outname+".log.txt","a")
             SIDTs.writelines("Begin of log file.\n")        
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.writelines("> Simulation took:{0:0.2f} seconds\n".format(reloj))
             SIDTs.writelines("> With # {0} reactions\n".format(str(numreacs)))
             SIDTs.writelines("> With mechanism : {0}\n".format(str(mecha2)))
             SIDTs.writelines("> Date of simulation : {0}\n".format(str(today)))
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.writelines("> System info : \n")
             SIDTs.writelines("> machine : {0}\n".format(str(machina)))
             SIDTs.writelines("> version : {0}\n".format(str(version)))
             SIDTs.writelines("> platform : {0}\n".format(str(platformA)))
             SIDTs.writelines("> uname : {0}\n".format(str(uname)))
             SIDTs.writelines("> system : {0}\n".format(str(system)))
             SIDTs.writelines("> # processors : {0}\n".format(str(procs)))
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.writelines("End of log file\n")
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.close()
        elif clock > 60:
             reloj = round((float(clock)/60),2)
             print("\t> Simulation took:{0:0.2f} minutes".format(reloj))
             print("Printing log.txt file...")
             SIDTs=open(outname+".log.txt","a")
             SIDTs.writelines("Begin of log file.\n")        
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.writelines("> Simulation took:{0:0.2f} minutes\n".format(reloj))
             SIDTs.writelines("> With # {0} reactions\n".format(str(numreacs)))
             SIDTs.writelines("> With mechanism : {0}\n".format(str(mecha2)))
             SIDTs.writelines("> Date of simulation : {0}\n".format(str(today)))
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.writelines("> System info : \n")
             SIDTs.writelines("> machine : {0}\n".format(str(machina)))
             SIDTs.writelines("> version : {0}\n".format(str(version)))
             SIDTs.writelines("> platform : {0}\n".format(str(platformA)))
             SIDTs.writelines("> uname : {0}\n".format(str(uname)))
             SIDTs.writelines("> system : {0}\n".format(str(system)))
             SIDTs.writelines("> # processors : {0}\n".format(str(procs)))
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.writelines("End of log file\n")
             SIDTs.writelines("--------------------------------------------\n")
             SIDTs.close()
        print("All done!.")
