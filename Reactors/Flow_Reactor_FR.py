###########################################
#  Created in 2018
#  @author: version 1.0. Sergio Martinez
###########################################
import os
import glob
import cantera as ct
import numpy as np
import re
import scipy
import scipy.optimize
import pandas as pd
import time
import datetime
import sys 
import time

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def InputFileReaderPFR(FileNamePath,OutputPath):
    FILENAME = FileNamePath.split('\\')
    INPUTFILENAME = FILENAME[-1]; 
    FilePFR= pd.read_csv(FileNamePath, names=['col'],sep='%s',engine='python');
    FolderName = FILENAME[3:-1]; 
    InputFileNames = (INPUTFILENAME.split(".PFR_EXP"))[0]; 
    FolderDir = "\\".join(FolderName) + "\\" + InputFileNames
    if not os.path.isdir(OutputPath + FolderDir):
           os.makedirs(OutputPath + FolderDir)
    FinalPath = str(OutputPath + FolderDir)
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # FILTERING DATA BY SPECIFIC NAMES, FROM HERE WE WILL GET THE FUEL,OXIDIZERS,
    # CONCENTRATIONS, INITIAL TEMPERATURES AND INITIAL PRESSURES, ETC... MAKING
    # DICTIONARIES TO STORAGE THE DATA READ IT.
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    NamesDict = {}
    ReactorsDict = {}
    UnitsDict = {}
    print('\t>',FileNamePath)
    try:
        REACS = (FilePFR[FilePFR.col.str.contains("REAC:")]).col.str.split(':')
    except:
        FilePFR = pd.read_fwf(FileNamePath,names=['col']);
        REACS  = (FilePFR[FilePFR.col.str.contains("REAC:")]).col.str.split(':')
    REACTANTS = []
    for i in range(len(REACS)):
        REAC1 = REACS.iloc[i][1]
        REAC1 = ":".join(REAC1.split())
        REACTANTS.append(REAC1)
    REACTANTS = str(",".join(REACTANTS))
    AUTHOR    = ((FilePFR[FilePFR.col.str.contains("AUTHOR:")]
               ).iloc[0].str.split(':'))[0][1]
    JOURNAL   = ((FilePFR[FilePFR.col.str.contains("JOURNAL:")]
                ).iloc[0].str.split(':'))[0][1]
    PDF       = ((FilePFR[FilePFR.col.str.contains("PDF:")]
            ).iloc[0].str.split(':'))[0][1]
    REACTOR   = ((FilePFR[FilePFR.col.str.contains("REACTOR:")]
                ).iloc[0].str.split(':'))[0][1]
    DATA      = ((FilePFR[FilePFR.col.str.contains("DATA:")]
             ).iloc[0].str.split(':'))[0][1]
    PRESSURE  = ((FilePFR[FilePFR.col.str.contains("PRESSURE:")]
                 ).iloc[0].str.split(':'))[0][1]
    #PRESSURE = str(PRESSURE)
    TIME      = ((FilePFR[FilePFR.col.str.contains("TIME:")]
            ).iloc[0].str.split(':'))[0][1]
    ENERGY_EQUATION = ((FilePFR[FilePFR.col.str.contains("ENERGY_EQUATION:")]
            ).iloc[0].str.split(':'))[0][1]
    #T_INITIAL = ((FilePFR[FilePFR.col.str.contains("T_INITIAL:")]
    #        ).iloc[0].str.split(':'))[0][1]
    #P_INITIAL = ((FilePFR[FilePFR.col.str.contains("P_INITIAL:")]
    #        ).iloc[0].str.split(':'))[0][1]
    #
    try:
        T_INITIAL = ((FilePFR[FilePFR.col.str.contains("T_INITIAL:")]
                          ).iloc[0].str.split(':'))[0][1]
        P_INITIAL = ((FilePFR[FilePFR.col.str.contains("P_INITIAL:")]
                          ).iloc[0].str.split(':'))[0][1]
        ReactorsDict["T_INITIAL"] = T_INITIAL
        ReactorsDict["P_INITIAL"] = P_INITIAL
    except:
        ReactorsDict["T_INITIAL"] = 0
        ReactorsDict["P_INITIAL"] = 0
    #
    PHI       = ((FilePFR[FilePFR.col.str.contains("PHI:")]
            ).iloc[0].str.split(':'))[0][1]
    TEMPERATURE = ((FilePFR[FilePFR.col.str.contains(
                "TEMPERATURE:")]).iloc[0].str.split(':'))[0][1]
    DropingStuff    = FilePFR[FilePFR.col.str.contains("!") == False]
    ITMM            = DropingStuff[DropingStuff.col.str.contains("TAU") == True]
    ITM             = pd.DataFrame((ITMM.col.str.split()).tolist())
    dataPFR         = pd.DataFrame((DropingStuff.col.str.split()).tolist())
    NamesDict["AUTHOR"]             = AUTHOR
    NamesDict["JOURNAL"]            = JOURNAL
    NamesDict["PDF"]                = PDF
    NamesDict["ENERGY_EQUATION"]    = ENERGY_EQUATION
    NamesDict["InpFileName"]        = INPUTFILENAME
    NamesDict["FolderDir"]          = FolderDir
    NamesDict["ITM"]                = ITM
    ReactorsDict["REACTANTS"]       = REACTANTS
    ReactorsDict["REACTOR"]         = REACTOR
    ReactorsDict["DATA"]            = DATA
    ReactorsDict["PHI"]             = PHI
    ReactorsDict["P_INITIAL"]       = P_INITIAL
    ReactorsDict["T_INITIAL"]       = T_INITIAL
    UnitsDict["PRESSURE"]           = PRESSURE
    UnitsDict["TEMPERATURE"]        = TEMPERATURE
    UnitsDict["TIME"]               = TIME
    return (NamesDict, ReactorsDict, UnitsDict,dataPFR,FinalPath,ITM)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def PFRCanteraSimulator(NamesDict,ReactorsDict,UnitsDict,dataPFR,FinalPath,ITM,MechFile):
    start     = time.time()
    TP        = dataPFR
    filename  = str(NamesDict["InpFileName"])
    NumCols   = len(TP.columns); print("\t> Number of species of interest: {0}".format(NumCols))

    LabelCols = TP.iloc[0,:]; INTER = [];
    for u in range(NumCols-1):
        INTER.append(str(ITM[u+1][0]))
        
    print("\t> List of species of interest:\n\t> {0}".format(INTER))
    ct.suppress_thermo_warnings(); now = datetime.datetime.now();
    NOW = str(now)[:10]; dash= '='*100; half= '-'*100; gato= '#'*100;
    #simname  = input('Please give a name for this job: ') or 'a'; 	 
    #print(simname)
    #########################
    # Input Parameters
    ########################
    T_0      = float(ReactorsDict["T_INITIAL"])
    # inlet temerature [K]
    pressure = float(ReactorsDict["P_INITIAL"])*ct.one_atm  # constant pressure [Pa]
    # ...................................................................
    if  T_0 == 0.0:
        Temps     = TP[0]
        Ti        = list(map(float, Temps.tolist()))
        Temps     = Ti 
        ResTime0  = TP[2]
        ResTimei  = [list(map(float, ResTime0))[-1]]
        resTime   = ResTimei
        press     = TP[1]
        Po        = press
        Pi        = [list(map(float, Po))[0]]
        pressures = Pi
        print(T_0, pressure)
    else:
        # .....................................................................
        print(T_0, pressure)
        # .......................................
        composition_0 = ReactorsDict["REACTANTS"]
        phi     = float(ReactorsDict["PHI"])
        length  = 8.0*5.0 # *approximate* PFR length [m]
        u_0     = 8.0  # inflow velocity [m/s]
        #area = 0.00812  # cross-sectional area [m**2]
        area    = 1.e-4  # cross-sectional area [m**2]
        # Resolution: 
        #The PFR will be simulated by 'n_steps' time steps or by a chain
        # of 'n_steps' stirred reactors.
        n_steps = 10000
        ######################## 
        ##########################################
        # Method 1: Lagrangian Particle Simulation
        ##########################################
        # A Lagrangian particle is considered which travels through the                     PFR. Its
        # state change is computed by upwind time stepping. The PFR result                       is produced
        # by transforming the temporal resolution into spatial locations.
        # The spatial discretization is therefore not provided a priori but                is instead
        # a result of the transformation.
        # import the gas model and set the initial conditions
        gas1 = ct.Solution(MechFile)
        mech2 = MechFile.rsplit('\\',1)[-1];
        mech3 = mech2.split(".cti")[0];
        print("\t>With {0} Mechanism\n ".format(mech2 )); print("="*60)
        gas1.TPX = T_0, pressure, composition_0
        mass_flow_rate1 = u_0 * gas1.density * area
        # create a new reactor
        r1 = ct.IdealGasConstPressureReactor(gas1)
        # create a reactor network for performing time integration
        sim1 = ct.ReactorNet([r1])
        # approximate a time step to achieve a similar resolution as in the next method
        t_total = 3100/T_0
        if float(TP.iloc[1,0]) == 0:
           dt = float(TP.iloc[2,0])*0.001
        else:
           dt = float(TP.iloc[1,0])*0.001 #t_total / n_steps
        # define time, space, and other information vectors
        #t1 = (np.arange(n_steps) + 1)* dt
        
        n1 = 0
        t  = 0.0
        #for n1, t_i in enumerate(t1):
        end_time = 20000 #int(((2*float(TP.iloc[-1,0])))/dt);
        t1 = np.zeros(end_time)
        z1 = np.zeros(end_time)
        u1 = np.zeros(end_time)
        states1 = ct.SolutionArray(r1.thermo)
        Tem=[]; 
        print("\t> This is your deltatime: {:.10f}".format(dt))
        print("\t> Simulation is going to finish when: {0} passed\n".format(end_time))
        print("-"*90)
        print("\t> Next file...")
        print("-"*90)
        while n1 < (end_time): # and r.T < T_equi - 0.1:
            t += dt
            # perform time integration
            sim1.advance(t)
            t1[n1] = sim1.time
            # compute velocity and transform into space
            u1[n1] = mass_flow_rate1 / area / r1.thermo.density
            z1[n1] = z1[n1 - 1] + u1[n1] * dt
            states1.append(r1.thermo.state)
            Tem.append(r1.T)
            n1 += 1
        for m in range(NumCols-1):
        
            if   NumCols == 2:
                 label1  = TP.iloc[0,1]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))]}); #sp1.append(df2)
            elif NumCols == 3:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))]}); #sp1.append(df2)
            elif NumCols == 4:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
                 label3  = TP.iloc[0,3]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))],str(label3):states1.X[:, gas1.species_index(str(label3))]}); #sp1.append(df2)
            elif NumCols == 5:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
                 label3  = TP.iloc[0,3]
                 label4  = TP.iloc[0,4]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))],str(label3):states1.X[:, gas1.species_index(str(label3))],str(label4):states1.X[:, gas1.species_index(str(label4))]}); #sp1.append(df2)
            elif NumCols == 6:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
                 label3  = TP.iloc[0,3]
                 label4  = TP.iloc[0,4]
                 label5  = TP.iloc[0,5]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))],str(label3):states1.X[:, gas1.species_index(str(label3))],str(label4):states1.X[:, gas1.species_index(str(label4))],str(label5):states1.X[:, gas1.species_index(str(label5))]}); #sp1.append(df2)
            elif NumCols == 7:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
                 label3  = TP.iloc[0,3]
                 label4  = TP.iloc[0,4]
                 label5  = TP.iloc[0,5]
                 label6  = TP.iloc[0,6]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))],str(label3):states1.X[:, gas1.species_index(str(label3))],str(label4):states1.X[:, gas1.species_index(str(label4))],str(label5):states1.X[:, gas1.species_index(str(label5))],str(label6):states1.X[:, gas1.species_index(str(label6))]}); #sp1.append(df2)
            elif NumCols == 8:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
                 label3  = TP.iloc[0,3]
                 label4  = TP.iloc[0,4]
                 label5  = TP.iloc[0,5]
                 label6  = TP.iloc[0,6]
                 label7  = TP.iloc[0,7]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))],str(label3):states1.X[:, gas1.species_index(str(label3))],str(label4):states1.X[:, gas1.species_index(str(label4))],str(label5):states1.X[:, gas1.species_index(str(label5))],str(label6):states1.X[:, gas1.species_index(str(label6))],str(label7):states1.X[:, gas1.species_index(str(label7))]}); #sp1.append(df2)
            elif NumCols == 9:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
                 label3  = TP.iloc[0,3]
                 label4  = TP.iloc[0,4]
                 label5  = TP.iloc[0,5]
                 label6  = TP.iloc[0,6]
                 label7  = TP.iloc[0,7]
                 label8  = TP.iloc[0,8]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))],str(label3):states1.X[:, gas1.species_index(str(label3))],str(label4):states1.X[:, gas1.species_index(str(label4))],str(label5):states1.X[:, gas1.species_index(str(label5))],str(label6):states1.X[:, gas1.species_index(str(label6))],str(label7):states1.X[:, gas1.species_index(str(label7))],str(label8):states1.X[:, gas1.species_index(str(label8))]}); #sp1.append(df2)
            elif NumCols == 10:
                 label1  = TP.iloc[0,1]
                 label2  = TP.iloc[0,2]
                 label3  = TP.iloc[0,3]
                 label4  = TP.iloc[0,4]
                 label5  = TP.iloc[0,5]
                 label6  = TP.iloc[0,6]
                 label7  = TP.iloc[0,7]
                 label8  = TP.iloc[0,8]
                 label9  = TP.iloc[0,9]
        
                 df2 = pd.DataFrame({str(label1):states1.X[:, gas1.species_index(str(label1))],str(label2):states1.X[:, gas1.species_index(str(label2))],str(label3):states1.X[:, gas1.species_index(str(label3))],str(label4):states1.X[:, gas1.species_index(str(label4))],str(label5):states1.X[:, gas1.species_index(str(label5))],str(label6):states1.X[:, gas1.species_index(str(label6))],str(label7):states1.X[:, gas1.species_index(str(label7))],str(label8):states1.X[:, gas1.species_index(str(label8))],str(label9):states1.X[:, gas1.species_index(str(label9))]}); #sp1.append(df2)
        SimTime = time.time() - start
        simname = NamesDict["InpFileName"]
        df1 = pd.DataFrame({'time':t1,'Ts':t1+0.01,'temperature':states1.T})
        df3 = pd.DataFrame({"REAC":ReactorsDict['REACTANTS'],
                            "PHI": ReactorsDict['PHI'],
                            "ORIGIN_FILE_NAME":simname,
                            "deltat /s":round(dt,7),
                            "SimTime /s":round((float(SimTime)),2), "DATE":NOW},index=[0])
        df  = df1.join(df2,lsuffix='_df1', rsuffix='_df2')
        dff = df.join(df3,lsuffix='_df', rsuffix='_df3')
        ds  = dff.to_csv(FinalPath+'\\'+mech3+'_PFR_EXP_.dat',index=False, sep="\t");
        # PLOTTING:
        ### delete the '#' symbol to activate the plotting features
        #for w in range(len(observables)):
        #    plt.rcParams["font.weight"] = "bold"
        #    plt.tick_params(axis= 'both', direction='out', length=4, width=2, labelsize=10)
        #    plt.rcParams["axes.labelweight"] = "bold"
        #    plt.rcParams["axes.linewidth"] = "3"
        #    #plt.scatter(ExpTempos[0],DFn[observables[w]], marker="s", label=str(observables[w].translate(subscript)), lw=4,color='black')
        #    plt.plot(ExpTempos[0],b[observables[w]], label=str(observables[w].translate(subscript)))
        #    #
        #    #
        #    #
        #    #plt.savefig(c + "\\" + NamesDict["InpFileName"] + "_" + observables[w] +'.png',dpi=800, bbox_inches="tight" )
        #    Y  = pd.to_numeric(DFn[observables[w]]) 
        #    YY = (Y).sort_values(ascending=False)
        #    plt.scatter(ExpTempos[0],YY, marker="s", label=str(observables[w].translate(subscript)), lw=4,color='black')
        #    plt.legend()
        #    plt.show()
        #    plt.close()

    return (ds)

#----------------------------------------------
#==============================================================================
if __name__ == "__main__":
   mech_list = [];
   path = "C:\\DKM\\MODELS\\"#os.getcwd();  
   fname = glob.glob(path+'/*.cti');
   #print('\t>Running with mechanism: \n\t',mechanism);
   #Mechanism = "/DKM/MODELS/19_48_C7.cti"
   for m in fname:
       print('\t>Running with mechanism: \n\t',m);
       Input = sys.argv[1] #"InputFileNameWithFullPath"
       (NamesDict, ReactorsDict, UnitsDict,dataPFR,FinalPath,ITM) = InputFileReaderPFR(Input, "C:\\DKM\\OUTPUTDATA\\");
       PFRCanteraSimulator(NamesDict, ReactorsDict, UnitsDict,dataPFR,FinalPath,ITM, m)
#=========================================================================
# Eppur si muove! - Galileo Galilei (1564-1642)
#=========================================================================
