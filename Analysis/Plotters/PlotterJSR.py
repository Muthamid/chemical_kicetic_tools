#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:23:58 2018
@author: Sergio Martinez
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
from matplotlib import colors as mcolors
from cycler import cycler
# -------------------------------------------------------------------
import statistics
from matplotlib.ticker import ScalarFormatter
# -------------------------------------------------------------------
def FileFinder(Expspath, ext):
    FILEFINDER  = []
    # traverse whole directory
    for root, dirs, files in os.walk(Expspath):
        # select file name
        for file in files:
            # check the extension of files
            if file.endswith(ext):
                # print whole path of files
                FILEFINDER.append(os.path.join(root, file))
            else:
                pass
    return FILEFINDER
# -------------------------------------------------------------------
# ---------------------------
counter  = 0
dicfiles = {}; FILES = [];
dicdirs  = {}; filesDirs = [];
for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            if ".JSR_EXP" in str(name):
                    FILES.append(os.path.join(root, name))
                    dicfiles[counter] = os.path.join(root, name)
                    counter += 1
# -----------------------------------------------------------------------------
print("\t> # of mechanism to plot: ",str(len(sys.argv)-1))
if    (len(sys.argv)) == 2:
      Mechas = [sys.argv[1]]
elif  (len(sys.argv)) == 3:
      Mechas = [sys.argv[1],sys.argv[2]]
elif  (len(sys.argv)) == 4:
      Mechas = [sys.argv[1],sys.argv[2],sys.argv[3]]
elif  (len(sys.argv)) == 5:
      Mechas = [sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]]
elif  (len(sys.argv)) == 6:
      Mechas = [sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]]
elif  (len(sys.argv)) == 7:
      print("\t> Error, plotter is only allowed to work with")
      print("\t> 3 mechanism names") 
      print("\t> Au revoir!")
      sys.exit(1)
else:
      print("\t> Error, plotter needs an argument...")
      print("\t> e.g. >python PloterRCM.py AramcoMech1.3") 
      print("\t> No extension cti/xml is required at the end of mechanism name")
      sys.exit(1)
# -----------------------------------------------------------------------------
print("\t> Mechanisms to plot", str(Mechas))
# -----------------------------------------------------------------------------
path = os.getcwd()
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
# Block of reading experiments:
# If you wanna  work with only one file, edit next line with the whole path of that file, and uncomment it (remove the dash symbol #)
#FILES = ["C:\\DKM\\DATA\\DATA_SHELL_EXP\\CH4\\SPECIES_TIME\\NATIVEL\\1_ATM_PYR_0.1.JSR_EXP"]
# If you wanna  work with only one folder, edit next line with the whole path of that directory, and uncomment it (remove the dash symbol #)
#FILES = glob.glob("C:\\DKM\\DATA\\DATA_SHELL_EXP\\CH4\\SPECIES_TIME\\NATIVEL\\*.JSR_EXP")
#FILES  = glob.glob("C:\\DKM\DATA\\DATA_SHELL_EXP\\CH4\\JSR\\NATIVEL\\*.JSR_EXP")
#print(FILES)
for i in FILES:
     FileJSR         = pd.read_csv(i,names=['col'],sep="%s",engine="python");
     DSpeciesJSR     = FileJSR[FileJSR.col.str.contains("!") == False]
     DSJSR           = DSpeciesJSR[DSpeciesJSR.col.str.contains("TAU") == True]
     DataEspeciesJSR = pd.DataFrame((DSJSR.col.str.split()).tolist())
     DropingStuff    = FileJSR[FileJSR.col.str.contains("!") == False]
     DropingStuff2   = DropingStuff[DropingStuff.col.str.contains("TAU") == False]
     dataJSR         = pd.DataFrame((DropingStuff2.col.str.split()).tolist())
     # **********************************************************************
     Species = list(DataEspeciesJSR)#.split('\\t')
     SP      = len(Species);
     Prods   = []; 
     for k in range(SP):
     #print(DataEspeciesJSR.iloc[0,k])
         if   str(DataEspeciesJSR.iloc[0,k]) == 'T/K':
              pass
         elif str(DataEspeciesJSR.iloc[0,k]) == 'P/ATM':
              pass
         elif str(DataEspeciesJSR.iloc[0,k]) == 'TAU/S':
              pass
         else:
              Prods.append(str(DataEspeciesJSR.iloc[0,k]))
     speciesToPlot = Prods
     #print(Species, SP, Prods, speciesToPlot, observables)
     # ***********************************************************************
     #### author and journal
     AUTHOR    = ((FileJSR[FileJSR.col.str.contains("AUTHOR:")]).iloc[0].str.split(':'))[0][1]
     JOURNAL   = ((FileJSR[FileJSR.col.str.contains("JOURNAL:")]).iloc[0].str.split(':'))[0][1]
     PDF       = ((FileJSR[FileJSR.col.str.contains("PDF:")]).iloc[0].str.split(':'))[0][1]			
     PHI       = ((FileJSR[FileJSR.col.str.contains("PHI:")]).iloc[0].str.split(':'))[0][1]			
     REACS     = (FileJSR[FileJSR.col.str.contains("REAC:")]).col.str.split(':')
     REACTANTS = []; OxDil = []; DILS = []; FuelNamesList = []; 
     print("\t Species: {0}".format(speciesToPlot))
     SUMAVAL = [];
     for m in range(len(REACS)):
            REAC1 = REACS.iloc[m][1]; 
            SUMAVAL.append(float(REAC1.split()[-1]))
     SumaRValues = sum(SUMAVAL); #print(SumaRValues)
     #*********************
     if SumaRValues < 2:
            for u in range(len(REACS)):
                if "C" in REACS.iloc[u][1]:
                   REAC1 = (REACS.iloc[u][1]); namb = round(float(REAC1.split()[-1])*100,3); FuelNamesList.append(REAC1.split()[0]);
                   REAC11 = str(namb)+"% "+str(REAC1.split()[-0]).translate(subscript) #":".join(REAC1.split())
                   REACTANTS.append(REAC11)
                if "O2" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1])*100,3);
                   REAC1 = str(namb)+"% "+str(REAC1.split()[-0]).translate(subscript) #":".join(REAC1.split())
                   OxDil.append(REAC1)
                if "N2" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1])*100,3);
                   if namb > 0:
                      REAC1 = str(namb)+"% "+str(REAC1.split()[-0]).translate(subscript) #":".join(REAC1.split())
                      DILS.append(REAC1)
                   else:
                      pass
                if "Ar" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1])*100,3);
                   if namb > 0:
                      REAC1 = str(namb)+"% Ar" #":".join(REAC1.split())
                      DILS.append(REAC1)
                   else:
                      pass
                if "AR" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1])*100,3);
                   if namb > 0:
                      REAC1 = str(namb)+"% Ar" #":".join(REAC1.split())
                      DILS.append(REAC1)
                if "CO2" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1])*100,3);
                   if namb > 0:
                      REAC1 = str(namb)+"% Ar" #":".join(REAC1.split())
                      DILS.append(REAC1)
                   else:
                      pass
     #*********************
     if SumaRValues >= 2:
            for u in range(len(REACS)):
                if "C" in REACS.iloc[u][1]:
                   REAC1 = (REACS.iloc[u][1]); namb = round(float(REAC1.split()[-1]),3); FuelNamesList.append(REAC1.split()[0]);
                   REAC11 = str(namb)+"% "+str(REAC1.split()[-0]).translate(subscript) #":".join(REAC1.split())
                   REACTANTS.append(REAC11)
                if "O2" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1]),3);
                   REAC1 = str(namb)+"% "+str(REAC1.split()[-0]).translate(subscript) #":".join(REAC1.split())
                   OxDil.append(REAC1)
                if "N2" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1]),3);
                   if namb > 0:
                      REAC1 = str(namb)+"% "+str(REAC1.split()[-0]).translate(subscript) #":".join(REAC1.split())
                      DILS.append(REAC1)
                   else:
                      pass
                if "Ar" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1]),3);
                   if namb > 0:
                      REAC1 = str(namb)+"% Ar" #":".join(REAC1.split())
                      DILS.append(REAC1)
                   else:
                      pass
                if "AR" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1]),3);
                   if namb > 0:
                      REAC1 = str(namb)+"% Ar" #":".join(REAC1.split())
                      DILS.append(REAC1)
                if "CO2" in REACS.iloc[u][1]:
                   REAC1 = REACS.iloc[u][1]; namb = round(float(REAC1.split()[-1]),3);
                   if namb > 0:
                      REAC1 = str(namb)+"% Ar" #":".join(REAC1.split())
                      DILS.append(REAC1)
                   else:
                      pass
     #*********************
     REACTANTS = str(", ".join(REACTANTS)); OxDil = str(", ".join(OxDil))
     DILS      = str(", ".join(DILS)); FNam = str("_".join(FuelNamesList))
     # ...........................................
     #print(REACTANTS, OxDil, DILS, FNam)
     print("\t> Plotting directory/file/case: ",i); #print(REACTANTS)
     temp     = pd.to_numeric(dataJSR.iloc[:,0]); 
     tempplot = pd.to_numeric(temp)
     PCC      = (dataJSR.iloc[:,1]); CompressedPressureAverage = [];
     for v in range(len(PCC)):
         CompressedPressureAverage.append(float(PCC.loc[v]))
     PC = round(statistics.mean(CompressedPressureAverage),2)
     for u in range(len(speciesToPlot)):
         # Activate latex text rendering
         Specie2plot1 = pd.to_numeric(dataJSR.iloc[:,u+3]); #print(tempplot,Specie2plot);
         plt.rcParams["font.weight"] = "bold"
         plt.tick_params(axis= 'both', direction='out', length=4, width=3, labelsize=10)
         plt.rcParams["axes.labelweight"] = "bold"
         plt.rcParams["axes.linewidth"] = "3"
         plt.rc('axes', prop_cycle=(cycler('color', ['k','r', 'b','g','m'])))
         legend_properties = {'weight':'bold'}
         plt.scatter(tempplot,Specie2plot1, marker="s", label=str(Prods[u].translate(subscript)), lw=4,color='black')
         FoldName  = ((i.rsplit("\\",1))[0]); #print(FoldName); print(Mechas) 
         FoldNamee = ((i.rsplit("\\",1))[-1]); FoldNamee = FoldNamee.split(".JSR_EXP")[0]
         LineStyles = ['-' , '--' , '-.' , ':'];
         # Reading Simulations:
         for j in Mechas:
              specs    = '\\'+str(FoldName); 
              FILESS   = glob.glob(FoldName+"\\"+j+"\\*\\"+str(FoldNamee)+".SP_OUT");
              if len(FILESS) > 0:
                 readmod        = pd.read_csv(FILESS[0],names=['col'],sep="%s",engine="python");
                 DropingStuff   = readmod[readmod.col.str.contains("!") == False]
                 dataJSRMod     = pd.DataFrame(DropingStuff.col.str.split().tolist())
                 dataJSRMod     = (dataJSRMod.rename(columns=dataJSRMod.iloc[0])).drop(dataJSRMod.index[0]) ; #dataST.columns = dataST.iloc[0]
                 legend_properties = {'weight':'bold'}
                 dataJSRMod     = dataJSRMod.reset_index(drop=True)
                 tempMod        = pd.to_numeric(dataJSRMod["T/K"]); #print(tempMod)
                 Specie2plot2   = pd.to_numeric((dataJSRMod.iloc[:,u+5]));#print(Specie2plot)
                 SpecieExpAverageValue = [];
                 for w in range(len(Specie2plot1)):
                     SpecieExpAverageValue.append(float(Specie2plot1.loc[w]))
                 SEAAV = statistics.mean(SpecieExpAverageValue)
                 CASE    = REACTANTS
                 # *********************
                 if len(REACTANTS) == 0:
                    title   = str(OxDil)+", "+str(DILS)+"\n$\phi$ = "+str(PHI.strip() )+", "+str(PC)+" atm";
                 if len(OxDil) == 0:
                    title   = str(CASE)+"\n "+str(DILS)+"\n$\phi$ = "+str(PHI.strip() )+", "+str(PC)+" atm";
                 if len(DILS) == 0:
                    title   = str(CASE)+"\n "+str(OxDil)+"\n$\phi$ = "+str(PHI.strip() )+", "+str(PC)+" atm";
                 if len(REACTANTS) > 0 and len(OxDil) > 0 and len(DILS) > 0:
                    title   = str(CASE)+"\n "+str(OxDil)+", "+str(DILS)+"\n$\phi$ = "+str(PHI.strip() )+", "+str(PC)+" atm";
                 # *********************
                 q = Mechas.index(j); #print(q); print(Mechas); print(j)
                 if     q ==0:
                      plt.plot(tempMod, Specie2plot2, linestyle= LineStyles[q], color="k", label=str(j), lw=3)
                 elif   q ==1:
                      plt.plot(tempMod, Specie2plot2, linestyle= LineStyles[q], color="k", label=str(j), lw=3)
                 elif   q ==2:
                      plt.plot(tempMod, Specie2plot2, linestyle= LineStyles[q], color="k", label=str(j), lw=3)
                 elif   q ==3:
                      plt.plot(tempMod, Specie2plot2, linestyle= LineStyles[q], color="k", label=str(j), lw=3)
                 plt.tick_params(axis = "y", which = "both", bottom = False, left = False)
                 plt.title(title, fontsize=14, weight= 'bold')
                 plt.xlabel('$\it{T}  /  K$', fontsize=20, weight= 'bold'); 
                 plt.ylabel('Mole fraction', fontsize=20, weight= 'bold')
                 plt.legend(fontsize="large",frameon=True)
                 plt.grid(b=True, which='major', color='grey', linestyle=':'); plt.grid(b=True, which='minor', color='grey', linestyle=':')
                 if SEAAV > 0:
                    MAX    = 1.2*max(Specie2plot1);
                    MIN    = 0.1*min(Specie2plot1);
                    plt.ylim(bottom=MIN,top=MAX)
                 elif SEAAV <= 0:
                    MAX    = 1.2*max(Specie2plot2);
                    MIN    = 0.1*min(Specie2plot2);
                    plt.ylim(bottom=MIN,top=MAX)
                 plt.tight_layout();
              elif len(FILESS) == 0:
                   pass
         directory = FoldName
         if not os.path.isdir(directory+"\\"+j+"\\JSRplots\\"):
                os.makedirs(directory+"\\"+j+"\\JSRplots\\")
         plt.savefig( directory+"\\"+j+"\\JSRplots\\"+"\\JSR_"+FNam+'_'+FoldNamee+"_"+str(Prods[u])+'.png',dpi=500)
         #plt.savefig( directory+"\\"+j+"\\JSRplots\\"+"\\"+FNam+'_'+FoldNamee+"_"+str(Prods[u])+'.png',dpi=500)
         plt.close()
# -------------------------------------------------------------------------------------------------
