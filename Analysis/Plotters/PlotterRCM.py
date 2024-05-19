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
counter  = 0
dicfiles = {}; FILES = [];
dicdirs  = {}; filesDirs = [];
for root, dirs, files in os.walk(".", topdown=False):
    for name in dirs:
        filesDirs.append(os.path.join(root, name))
        dicdirs[counter] = os.path.join(root, name)
        counter += 1
for root, dirs, files in os.walk(".", topdown=False):
    for k in filesDirs:
        File2track = k.rsplit('\\',1)[-1]; #print(File2track)
        for name in files:
            #print(os.path.join(root, name))
            if "NUIGMech1.1" in str(name):
                 pass
            elif "21_08_C4" in str(name):
                 pass
            elif "21_08_C4_ahmed" in str(name):
                 pass
            elif File2track+".RCM_IDT" == str(name):
            #if File2track+".png" == str(name):
                FILES.append(os.path.join(root, name))
                dicfiles[counter] = os.path.join(root, name)
                counter += 1
#print(os.path.join(root, name))
#print(dirs)
#print("LOL")
#print("FILES:")
#print(FILES)
#print("DIRS:")
#print(filesDirs)
#Mechas = ["20_11_C4","20_11_C4_100atmRate"]
#print(len(sys.argv))
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

#------------------------------------------------------------------------------
print("\t> Mechanisms to plot", str(Mechas))
path = os.getcwd()
# EXPERIMENTS
#FILES = glob.glob(path+'\\EXP.dat');
#for p in range(2):
# READING EXPERIMENTAL INPUT FILE
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
for i in FILES:
    readexp       = pd.read_csv(i,names=['col'],sep="%s",engine="python");
    DropingStuff  = readexp[readexp.col.str.contains("!") == False]
    DropingStuff2 = DropingStuff[DropingStuff.col.str.contains("TC") == False]
    dataRCM       = pd.DataFrame(DropingStuff2.col.str.split().tolist())
    #### author and journal
    AUTHOR  = ((readexp[readexp.col.str.contains("AUTHOR:")]).iloc[0].str.split(':'))[0][1]
    JOURNAL = ((readexp[readexp.col.str.contains("JOURNAL:")]).iloc[0].str.split(':'))[0][1]
    PDF = ((readexp[readexp.col.str.contains("PDF:")]).iloc[0].str.split(':'))[0][1]			
    REACS  = (readexp[readexp.col.str.contains("REAC:")]).col.str.split(':')
    # ----------------------------------------------------------------------------------------------------
    REACTANTS     = []
    OxDil         = []
    DILS          = []
    FuelNamesList = []
    # --------------------------------------------------------------------
    SUMAVAL = []; #print(TimeScale);
    for m in range(len(REACS)):
           REAC1 = REACS.iloc[m][1]; 
           SUMAVAL.append(float(REAC1.split()[-1]))
    SumaRValues = sum(SUMAVAL); #print(SumaRValues)
    # ---------------------------------------------------------------------
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
    # ----------------------------
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
    # ----------------------------
    REACTANTS = str(", ".join(REACTANTS)); OxDil = str(", ".join(OxDil))
    DILS      = str(", ".join(DILS)); 
    ####################################################
    #for u in range(len(REACS)):
    #    if "N2" in REACS[u]:
    #        print(str(REACS[u]))
    ####
    print("\t> Plotting directory/file/case: ",i)
    temp     = pd.to_numeric(dataRCM.iloc[:,2]); 
    idt      = pd.to_numeric(dataRCM.iloc[:,5]);
    tempplot = 1000/temp; 
    idtplot  = (idt); Yerr = [(0.2*x) for x in idt] 
# activate latex text rendering 
    plt.rcParams["font.weight"] = "bold"
    plt.tick_params(axis= 'both', direction='out', length=4, width=3, labelsize=10)
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams["axes.linewidth"] = "3"
    plt.rc('axes', prop_cycle=(cycler('color', ['k','r', 'b','g','m'])))
    legend_properties = {'weight':'bold'}
    plt.yscale("log")
    plt.scatter(tempplot,idtplot, marker="s", label='Exp Data', lw=4,color='black')
    plt.errorbar(tempplot,idtplot, yerr=Yerr, label=None, fmt='none', ecolor="k" )
# SIMULATIONS
    for j in Mechas:
        FoldName = ((i.rsplit('\\',1))[0]); FNam = ((i.split('\\'))[-5]); #print(FoldName)
        specs    = '\\'+str(FoldName); 
        FILESS   = glob.glob(FoldName+'\\'+j+'.RCM_IDT');
        #print(FILESS)
        plt.yscale("log")
        for m in range(len(FILESS)):
            legend_properties = {'weight':'bold'}
            readmod   = pd.read_table(FILESS[m]);
            mech      = (((FILESS[m].split('\\'))[-1]).split('.RCM_IDT'))[0];
            tempMod2  = (readmod.iloc[:,4]); 
            tempMod2  = [x for x in tempMod2]; tempMod2 = sorted(tempMod2,reverse=True)
            idtTMod2  = (readmod.iloc[:,6]); #readmod['IDT/s'],reverse=True);
            idtTMod2  = [1000*(x) for x in idtTMod2]
            idtTMod2  = sorted((idtTMod2),reverse=True)
            idtMod2   = (readmod.iloc[:,5]); #readmod['IDT/s'],reverse=True);
            if readmod.iloc[0,5] > 1:
               idtMod2   = [(x) for x in idtMod2]
            else:
               idtMod2   = [1000*(x) for x in idtMod2]
            idtMod2   = sorted((idtMod2),reverse=True)
            #FUEL    = readmod.iloc[0,9]; FUEL1st = (str(FUEL).split(":"))[0]
            #print(tempMod2,idtMod2); 
            tempMod = []; idtMod = [];
            for m in range(len(idtMod2)):
                if idtMod2[m] > 0:
                   idtMod.append(idtMod2[m])
                   tempMod.append(tempMod2[m])
                elif idtMod2[m] == 0 and idtTMod2[m] > 0:
                   idtMod.append(idtTMod2[m])
                   tempMod.append(tempMod2[m])
                elif idtMod2[m] < 0:
                   idtMod.append(idtTMod2[m])
                   tempMod.append(tempMod2[m])
                elif idtMod2[m] < 0 and idtTMod2[m] < 0:
                   idtMod.append(0.0)
                   tempMod.append(tempMod2[m])

            CASE    = REACTANTS # i.split("\\")[1];    
            initconds = FoldName.split("\\")[-1]; PHIs = initconds.split("_")[0]; 
            PRESss  = (initconds.split("_")[1]).rsplit("ATM",1)[0]; PRESs = (PRESss)      
            if   PHIs == "F050":
                 title   = str(CASE)+"\n "+str(OxDil)+", "+str(DILS)+"\n $\phi$ = 0.5, "+str(PRESs)+" atm";
            elif PHIs == "F100":
                 title   = str(CASE)+"\n"+str(OxDil)+", "+str(DILS)+"\n $\phi$ = 1.0, "+str(PRESs)+" atm";
            elif PHIs == "F200":
                 title   = str(CASE)+"\n"+str(OxDil)+", "+str(DILS)+"\n $\phi$ = 2.0, "+str(PRESs)+" atm";
            else:
                 title   = str(CASE)+"\n"+str(OxDil)+", "+str(DILS)+"\n" +str(PHIs)+", "+str(PRESs)+" atm";
            plt.yscale("log");
            if   str(mech) == str(Mechas[0]):
                 print("Yes, {0} mechanism found!".format(mech)); mecham = mech;
                 mecham      = mech
                 if   mecham == "N13_C4_gas":
                      mecham = "NUIGMech1.3"
                 elif mecham == "NUIGMech1.2_C4":
                      mecham = "NUIGMech1.2"
                 elif mecham == "AM3.0":
                      mecham = "AramcoMech3.0"
                 else:               
                      pass
                 print(mecham)
                 mech = mecham
                 plt.plot(tempMod, idtMod, linestyle= '-', color="k", label=str(mech), lw=3)            
            elif str(mech) == str(Mechas[1]):
                 print("Yes, NUIGMech1.1 mechanism found!");
                 mecham      = mech
                 if   mecham == "N13_C4_gas":
                      mecham = "NUIGMech1.3"
                 elif mecham == "NUIGMech1.2_C4":
                      mecham = "NUIGMech1.2"
                 elif mecham == "AM3.0":
                      mecham = "AramcoMech3.0"
                 else:               
                      pass
                 print(mecham)
                 mech = mecham
                 plt.plot(tempMod, idtMod, linestyle= '-',color="r", label=str(mech), lw=3)
            else:
                 mecham      = mech
                 if   mecham == "N13_C4_gas":
                      mecham = "NUIGMech1.3"
                 elif mecham == "NUIGMech1.2_C4":
                      mecham = "NUIGMech1.2"
                 elif mecham == "AM3.0":
                      mecham = "AramcoMech3.0"
                 else:               
                      pass
                 print(mecham)
                 mech = mecham
                 plt.plot(tempMod, idtMod, linestyle= '--',label=str(mech), lw=3)
            plt.title(title, fontsize=14, weight= 'bold')
            plt.xlabel('1000 K / $\it{T}$', fontsize=20, weight= 'bold'); 
            plt.ylabel('IDT / ms', fontsize=20, weight= 'bold')
            #plt.ylim( bottom = 0.5*min(idtplot),top=2*max(idtplot)) #bottom
            plt.legend(frameon=True)#frameon=False); #plt.tight_layout()
    plt.grid(b=True, which='major', color='grey', linestyle=':'); plt.grid(b=True, which='minor', color='grey', linestyle=':')
    #if   "BAIGMOHOMMADI" in AUTHOR:
    #     plt.annotate(str(AUTHOR)+" - NUIG", (0.5,0.3), (0, -50),bbox=dict(facecolor='white', alpha=0.8), xycoords='axes fraction', textcoords='offset points', va='top',)
    #elif "RAMALINGAM" in AUTHOR:
    #     plt.annotate(str(AUTHOR)+" - AACHEN", (0.5,0.3), (0, -50),bbox=dict(facecolor='white', alpha=0.8), xycoords='axes fraction', textcoords='offset points', va='top')
    #else:
    #     plt.annotate(str(AUTHOR), (0.5,0.3), (0, -50),bbox=dict(facecolor='white', alpha=0.8), xycoords='axes fraction', textcoords='offset points', va='top')
    #plt.annotate(str(AUTHOR), (0.5,0.3), (0, -50), xycoords='axes fraction', textcoords='offset points', va='top')
    #plt.annotate(str(PDF), (0.6,0.4), (0, -50), xycoords='axes fraction', textcoords='offset points', va='top')
    #plt.text(0.98*MaxX, MaxY, 'Solid line: '+str(mecham),horizontalalignment='left',verticalalignment='top', fontsize=14)
    directory = FoldName.rsplit("\\",1)[0]; FoldNamee = FoldName.rsplit("\\",1)[-1]
    plt.tight_layout()
    plt.savefig( directory+"\\"+FNam+'_'+FoldNamee+'.png',dpi=500)
    plt.close()
    #plt.show()
#
path = os.getcwd()
# Get List of all images
FILES = glob.glob(path + '\*\RCM\*\*.png', recursive=True)
# ---
# Block for monofuels:
# ---
###f = open("RCM_data.LateX.out", "w")
####print((FILES))
###for k in FILES:
###    FNam         = (k.split('\\'))[-4]; #print(FNam+"\n") 
###    FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
###    EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
###    #print(k)
###    kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
###    f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
###f.close()
MONOS        = ["CH4","C2H4","C2H6","C3H8"]
BINARIES     = ["CH4_C2H4","CH4_C2H6","C2H4_C2H6","C2H4_C3H8","C2H6_C3H8"]
TERNARIES    = ["CH4_C2H4_C2H6"]
QUATERNARIES = ["CH4_C2H4_C2H6_C3H8"]
#print(FILES[0])
# ******
# MONOS
# ******
f = open("RCM_data.LateX.out", "w")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == MONOS[0]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == MONOS[1]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == MONOS[2]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == MONOS[3]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
# -------------------------------------------------------------------------------------------------
# ******
# BINARIES
# ******
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == BINARIES[0]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == BINARIES[1]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == BINARIES[2]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == BINARIES[3]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == BINARIES[4]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
# -------------------------------------------------------------------------------------------------
# ******
# TERNARIES
# ******
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == TERNARIES[0]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
# -------------------------------------------------------------------------------------------------
# ******
# QUATERNARIES
# ******
f = open("RCM_data.LateX.out", "a")
for k in FILES:
    FNam         = ((k.split('\\'))[-4]); #print(FNam)
    if    str(FNam) == QUATERNARIES[0]:
         FORCE        = (k.rsplit('\\',1))[0]; #print(FORCE+"\n")
         EachSplitted = ((k.split('\\'))[-1]).split(".png")[0]; #print(EachSplitted+"\n") 
         kiss = (FORCE+"/{"+EachSplitted+"}.png"); newPath = kiss.replace(os.sep, '/')
         f.write("\includegraphics[scale=0.5]{"+str(newPath)+"}\n")
f.close()
