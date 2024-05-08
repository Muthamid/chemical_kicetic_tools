# -*- coding: utf-8 -*-
#"""
#Created on 2018
# #######################################################
# # Sergio Martinez C3-NUIG, IE. 12-12-2019            ##
# # S1E9S8M1J1M9S9P4, "La verdad os hara libres!"      ##
# #######################################################
# This code snippet runs a LBV simulation using Cantera.
# """
import os
import sys
import glob
import cantera as ct
import numpy as np
from matplotlib import pyplot as plt
import re
import scipy
import scipy.optimize
import pandas as pd
import time
import datetime
# ##############
# 
# ##############

def InputFileReaderLBV(FileNamePath,OutputPath):
	FILENAME = FileNamePath.split('\\')
	INPUTFILENAME = FILENAME[-1]; 
	FileLBV= pd.read_csv(FileNamePath, names=['col'],sep='%s',engine='python');
	FolderName = FILENAME[3:-1]; 
	InputFileNames = (INPUTFILENAME.split(".PFP_EXP"))[0]; 
	FolderDir = "\\".join(FolderName) + "\\" + InputFileNames
	if not os.path.isdir(OutputPath + FolderDir):
		   os.makedirs(OutputPath + FolderDir)
	FinalPath = str(OutputPath + FolderDir)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FILTERING DATA BY SPECIFIC NAMES, FROM HERE WE WILL GET THE FUEL,OXIDIZERS,
# CONCENTRATIONS, INITIAL TEMPERATURES AND INITIAL PRESSURES, ETC... MAKING
# DICTIONARIES TO STORAGE THE DATA READ IT.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	NamesDict = {}
	ReactorsDict = {}
	UnitsDict = {}
	print('\t>',FileNamePath)
	AUTHOR = ((FileLBV[FileLBV.col.str.contains("AUTHOR:")]
               ).iloc[0].str.split(':'))[0][1]
	JOURNAL = ((FileLBV[FileLBV.col.str.contains("JOURNAL:")]
                ).iloc[0].str.split(':'))[0][1]
	PDF = ((FileLBV[FileLBV.col.str.contains("PDF:")]
            ).iloc[0].str.split(':'))[0][1]
	REACTOR = ((FileLBV[FileLBV.col.str.contains("REACTOR:")]
                ).iloc[0].str.split(':'))[0][1]
	DATA = ((FileLBV[FileLBV.col.str.contains("DATA:")]
             ).iloc[0].str.split(':'))[0][1]
	PRESSURE = ((FileLBV[FileLBV.col.str.contains("PRESSURE:")]
                 ).iloc[0].str.split(':'))[0][1]
	#PRESSURE = str(PRESSURE)
	#COLLECT = ((FileLBV[FileLBV.col.str.contains("COLLECTED_BY:")]
    #        ).iloc[0].str.split(':'))[0][1]
	#DATE = ((FileLBV[FileLBV.col.str.contains("DATE:")]
    #        ).iloc[0].str.split(':'))[0][1]
	#SOURCE = ((FileLBV[FileLBV.col.str.contains("SOURCE:")]
    #        ).iloc[0].str.split(':'))[0][1]
	#INTM = ((FileLBV[FileLBV.col.str.contains("INTM:")]
    #        ).iloc[0].str.split(':'))[0][1]
	TEMPERATURE = ((FileLBV[FileLBV.col.str.contains("TEMPERATURE:")]).iloc[0].str.split(':'))[0][1]
	DropingStuff    = FileLBV[FileLBV.col.str.contains("!") == False]
	TimeScale       = ((FileLBV[FileLBV.col.str.contains("LBV:")]).iloc[0].str.split(':'))[0][1]
	dataLBV         = pd.DataFrame((DropingStuff.col.str.split()).tolist())
	NamesDict["AUTHOR"]             = AUTHOR
	NamesDict["JOURNAL"]            = JOURNAL
	NamesDict["PDF"]                = PDF
	#NamesDict["COLLECT"]            = COLLECT
	#NamesDict["DATE"]               = DATE
	#NamesDict["SOURCE"]             = SOURCE
	NamesDict["InpFileName"]        = INPUTFILENAME
	NamesDict["FolderDir"]          = FolderDir
	ReactorsDict["REACTOR"]         = REACTOR
	ReactorsDict["DATA"]            = DATA
	#ReactorsDict["INTM"]            = INTM
	UnitsDict["PRESSURE"]           = PRESSURE
	UnitsDict["TEMPERATURE"]        = TEMPERATURE
	UnitsDict["TIMESCALE"]          = TimeScale
	return (NamesDict, ReactorsDict, UnitsDict,dataLBV,FinalPath)
# ####################################

OUTPUT_PATH= "C:\\DKM\\OUTPUTDATA\\"
dash= '='*90; half= '-'*90; gato = '#'*90; 

# ####################################
def LBVCanteraSimulator(MechFile, Input):
    start = time.time()
    (NamesDict, ReactorsDict, UnitsDict, TP, FinalPath) = InputFileReaderLBV(Input,OUTPUT_PATH)
    filename = str(NamesDict["InpFileName"])
    ct.suppress_thermo_warnings(); now = datetime.datetime.now();
    NOW = str(now)[:10]; #INTM = ((ReactorsDict["INTM"]).split(' '));
    #print(TP)
    #INTM2 = []; 
    #for p in range(len(INTM)):
    #    if INTM[p] == " ":
    #       continue
    #    else:
    #       INTM2.append(INTM[p])
    if not os.path.isdir(FinalPath + '\\SENS\\'):
       os.makedirs(FinalPath + '\\SENS\\'); 
    #assert os.path.exists(cantera_file_path)
    gas = ct.Solution(MechFile); #print(INTM);#print(INTM2);
    mech2 = MechFile.rsplit('\\',1)[-1];
    mech3 = mech2.split(".cti")[0];
    Ti = []; Pi = []; PHI = []; LVB = []; SPEED = [];
    ERROR = []; FITTED = [];
    print("\t>With {0} Mechanism\n ".format(mech2 )); print("="*60)
    CHMols = []; O2Mols = []; DilMols = []; ProdMols = [];
    for i in (TP.iloc[0,:]):
        if "CO2" in i or "H2O" in i:
             ProdMols.append(i)
        elif "N2" in i or "AR" in i or "HE" in i:
             DilMols.append(i)
        elif "O2" in i:
             O2Mols.append(i)
        elif "C" in i :
             CHMols.append(i)
        elif "H2" in i :
             CHMols.append(i)
    if len(CHMols) == 1:
       fuel = "mono"
    elif len(CHMols) == 2:
       fuel = "binary"
    elif len(CHMols) == 6:
       fuel = "puto"
    else:
       print("More than two and puto fuels")
       sys.exit(1)
    print(len(CHMols),CHMols)
    NumCols   = len(TP.columns); #print(NumCols)
    LabelCols = TP.iloc[0,:]; #print(LabelCols[6])
    for f in range(len(TP)-1):
        if fuel == "mono":
           fuels    = CHMols[0]+':'+str(TP.iloc[(f+1),4])
           oxidizer = O2Mols[0]+':'+str(TP.iloc[(f+1),5])
           dilution = DilMols[0]+':'+str(TP.iloc[(f+1),6])
        elif fuel == "binary":
           fuels    = CHMols[0]+':'+str(TP.iloc[(f+1),4])+","+CHMols[1]+':'+str(TP.iloc[(f+1),5])
           oxidizer = O2Mols[0]+':'+str(TP.iloc[(f+1),6])
           dilution = DilMols[0]+':'+str(TP.iloc[(f+1),7])
        elif fuel == "puto":
           fuels    = CHMols[0]+':'+str(TP.iloc[(f+1),4])+","+CHMols[1]+':'+str(TP.iloc[(f+1),5])+","+CHMols[2]+':'+str(TP.iloc[(f+1),6])+","+CHMols[3]+':'+str(TP.iloc[(f+1),7])+","+CHMols[4]+':'+str(TP.iloc[(f+1),8])+","+CHMols[5]+':'+str(TP.iloc[(f+1),9])
           oxidizer = O2Mols[0]+':'+str(TP.iloc[(f+1),6])
           dilution = DilMols[0]+':'+str(TP.iloc[(f+1),7])

        resto       = oxidizer+','+dilution
        ExpResult   = str(TP.iloc[(f+1),3]); LVB.append(ExpResult)
        phi         = float(TP.iloc[(f+1),2]); PHI.append(str(phi))
        temperature = float(TP.iloc[(f+1),1]); Ti.append(str(temperature))
        pressure    = float(TP.iloc[(f+1),0]); Pi.append(str(pressure))
        print("T = {} K".format(temperature))
        print("P = {} atm".format(pressure))
        gas.set_equivalence_ratio(phi, fuels, resto)
        gas.mole_fraction_dict()
        gas.TP = (temperature,pressure*ct.one_atm)
        gas()
        print(gas())
        #++++++++++++++++++++++++++
        #
        #++++++++++++++++++++++++++
        def extrapolate_uncertainty(grids, speeds,f):
            """
            Given a list of grid sizes and a corresponding list of flame speeds,
            extrapolate and estimate the uncertainty in the final flame speed.
            Also makes a plot.
            """
            grids = list(grids)
            speeds = list(speeds)
            #+++++++++++++++++++++++
            #
            #+++++++++++++++++++++++
	
            def speed_from_grid_size(grid_size, true_speed, error):
                """
                Given a grid size (or an array or list of grid sizes)
                return a prediction (or array of predictions)
                of the computed flame speed, based on 
                the parameters `true_speed` and `error`
                """
                return true_speed +  error * np.array(grid_size)**-1.
	
            popt, pcov = scipy.optimize.curve_fit(speed_from_grid_size, grids[-4:], speeds[-4:])
            perr = np.sqrt(np.diag(pcov))
            true_speed  = popt[0]
            percent_error_in_true_speed = 100.0*perr[0] / popt[0]
            print("Fitted true_speed is {:.4f} ± {:.4f} cm/s ({:.1f}%)".format(popt[0]*100,perr[0]*100,percent_error_in_true_speed))
            Fitted_speed = str(round(popt[0]*100,3))+'±'+str(round(perr[0]*100,3))+'('+str(round(percent_error_in_true_speed,3))+')'
            #print "convergerce rate wrt grid size is {:.1f} ± {:.1f}".format(popt[2], perr[2])
            estimated_percent_error = 100. * (speed_from_grid_size(grids[-1], *popt) - true_speed) / true_speed
            print("Estimated error in final calculation {:.1f}%".format(estimated_percent_error))
            total_error_estimate = abs(percent_error_in_true_speed) + abs(estimated_percent_error)
            print("Estimated total error {:.1f}%".format(total_error_estimate))
            #+++++++++++++++++++++++++++++++++++++++++++++++
            # Plotting
            #+++++++++++++++++++++++++++++++++++++++++++++++
            ###plt.semilogx(grids,speeds,'o-')
            ###plt.ylim(min(speeds[-5:]+[true_speed-perr[0]])*.95, max(speeds[-5:]+[true_speed+perr[0]])*1.05)
            ###plt.plot(grids[-4:], speeds[-4:], 'or')
            ###extrapolated_grids = grids + [grids[-1] * i for i in range(2,8)]
            ###plt.plot(extrapolated_grids,speed_from_grid_size(extrapolated_grids,*popt),':r')
            ###plt.xlim(*plt.xlim())
            ###plt.hlines(true_speed, *plt.xlim(), colors=u'r', linestyles=u'dashed')
            ###plt.hlines(true_speed+perr[0], *plt.xlim(), colors=u'r', linestyles=u'dashed', alpha=0.3)
            ###plt.hlines(true_speed-perr[0], *plt.xlim(), colors=u'r', linestyles=u'dashed', alpha=0.3)
            ###plt.fill_between(plt.xlim(), true_speed-perr[0],true_speed+perr[0], facecolor='red', alpha=0.1 )
            ####plt.text(grids[-1],speeds[-1],"{:.1f}%".format(estimated_percent_error))
            ###above = popt[1]/abs(popt[1]) # will be +1 if approach from above or -1 if approach from below
            ###plt.annotate("",
            ###        xy=(grids[-1], true_speed),
            ###        xycoords='data',
            ###        xytext=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
            ###         textcoords='data',
            ###         arrowprops=dict(arrowstyle='|-|',
            ###                    connectionstyle='arc3',
            ###                    color='black', shrinkA=0, shrinkB=0),
            ###        )
            ###plt.annotate("{:.2f}%".format(abs(estimated_percent_error)),
            ###    xy=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
            ###     xycoords='data',
            ###    xytext=(10,20*above),
            ###     textcoords='offset points',
            ###     arrowprops=dict(arrowstyle='->',
            ###                    connectionstyle='arc3')
            ###    )
            ###plt.annotate("",
            ###    xy=(grids[-1]*4, true_speed-(above*perr[0])),
            ###     xycoords='data',
            ###    xytext=(grids[-1]*4, true_speed),
            ###     textcoords='data',
            ###     arrowprops=dict(arrowstyle='|-|',
            ###                    connectionstyle='arc3',
            ###                    color='black', shrinkA=0, shrinkB=0),
            ###    )
            ###plt.annotate("{:.2f}%".format(abs(percent_error_in_true_speed)),
            ###    xy=(grids[-1]*4, true_speed-(above*perr[0])),
            ###     xycoords='data',
            ###    xytext=(10,-20*above),
            ###     textcoords='offset points',
            ###     arrowprops=dict(arrowstyle='->',
            ###                    connectionstyle='arc3')
            ###    )
            ###plt.ylabel("Flame speed (m/s)")
            ###plt.xlabel("Grid size")
            ###plt.savefig(FinalPath+'\\IMAGES\\'+filename+'_PHI_'+str(phi)+
            ###            '_'+str(mech3)+'_'+str(f)+'_extrapolate.png',dpi=800)
            ###plt.close()
            return (true_speed, total_error_estimate,Fitted_speed)
	
        #++++++++++++++++
        #
        #++++++++++++++++
        def make_callback(flame):
            speeds = []
            grids  = []
            #+++++++++++
            #
            #+++++++++++
            def callback(_):
                speed = flame.u[0]
                grid = len(flame.grid)
                speeds.append(speed)
                grids.append(grid)
                print("Iteration {}".format(len(grids)))
                print("Current flame speed is is {:.4f} cm/s".format(speed*100.))
                if len(grids) < 5:
                   #+++++++++++++++++++++++++++++++++++++++++++++++
                   # Returnning data
                   #+++++++++++++++++++++++++++++++++++++++++++++++
                    return 1.0 # 
                try:
                    extrapolate_uncertainty(grids, speeds,f)
                except Exception as e:
                    print("Couldn't estimate uncertainty", e.message)
                   #+++++++++++++++++++++++++++++++++++++++++++++++
                   # Returnning data
                   #+++++++++++++++++++++++++++++++++++++++++++++++
                    return 1.0 # continue anyway
                #+++++++++++++++++++++++++++++++++++++++++++++++
                # Returnning data
                #+++++++++++++++++++++++++++++++++++++++++++++++
                return 1.0
            #+++++++++++++++++++++++++++++++++++++++++++++++
            # Returnning data
            #+++++++++++++++++++++++++++++++++++++++++++++++
            return callback, speeds, grids
        # flame.set_steady_callback(make_callback()[0])
        # Domain width in metres
        width = 0.015
        # Create the flame object
        flame = ct.FreeFlame(gas, width=width)
        # Define tolerances for the solver
        # (these are used throughout the notebook)
        #refine_criteria = {'ratio':3, 'slope': 0.1, 'curve': 0.1} # around 145 grid size, Error from 0.4 - 6% ## refine_criteria = {'ratio':2, 'slope': 0.01, 'curve': 0.01} # around 1000 grid size, Error from 0.0 - 0.6% 
        refine_criteria = {'ratio':2, 'slope': 0.01, 'curve': 0.01}
        flame.set_refine_criteria(**refine_criteria)
        #flame.set_max_grid_points(flame.domains[flame.domain_index('flame')], 1e4)
        try:
           flame.set_max_grid_points(flame.domains[flame.domain_index('flame')], 1e4)
        except:
           flame.set_max_grid_points(flame.domains[flame.domain_index('flame')], 1e6)
        callback, speeds, grids = make_callback(flame)
        flame.set_steady_callback(callback)
        # Define logging level
        loglevel = 1
        flame.solve(loglevel=loglevel, auto=True)
        final_true_speed, percentage_uncertainty,Fitted_speed = extrapolate_uncertainty(grids, speeds,f)
        print("Final grid was size {}".format(grids[-1]))
        print("Final speed was {:.4f} cm/s".format(100*speeds[-1]))
        print("Estimated uncertainty is {:.1f}%".format(percentage_uncertainty))
        print("i.e. {:.3f} +/- {:.3f} cm/s".format(100*speeds[-1],
                                           percentage_uncertainty*speeds[-1]))
        for i in range(4,len(grids)):
            print("At step {}".format(i))
            print("Grid was size {}".format(grids[i]))
            print("Speed was {:.4f} cm/s".format(100*speeds[i]))
            true_speed, percentage_uncertainty,Fitted_speed = extrapolate_uncertainty(grids[:i], speeds[:i],f)
            print("Estimated uncertainty was {:.1f}%".format(percentage_uncertainty))
            print("i.e. {:.3f} +/- {:.3f} cm/s".format(100*speeds[i],
                                           percentage_uncertainty*speeds[i]))
            print("or  {:.3f} -- {:.3f} cm/s".format((100-percentage_uncertainty)*speeds[i],
                                           (100+percentage_uncertainty)*speeds[i]))
            print("(For reference, the 'final' extrapolated speed was {:.3f} cm/s".format(100*final_true_speed))
            print("="*80)
        SPEED.append(round(100*final_true_speed,3))
        ERROR.append(round(percentage_uncertainty,3))
        FITTED.append(Fitted_speed)
        ## Solve with multi-component transport properties
        flame.transport_model = 'Multi'
        flame.soret_enabled = True
        flame.solve(loglevel) # don't use 'auto' on subsequent solves
        ##print('multicomponent flamespeed = {0} m/s'.format(flame.u[0]))
    ####+++++++++++++++++++++++++++++++++++++++++++++++
    #### Plotting sensitivity analyses
    ####+++++++++++++++++++++++++++++++++++++++++++++++
    if "sens" in sys.argv:
        #+++++++++++++++++++++++++++++++++++++++++++++++
        # Starting sensitivity analysis of the reactions
        #+++++++++++++++++++++++++++++++++++++++++++++++
        #ListOfEqsRats = [0.5, 1.0, 2.0]
        sens1 = flame.get_flame_speed_reaction_sensitivities()
        # print the reaction number followed by the associated sensitivity 
        #print()
        #print('Rxn #   k/S*dS/dk    Reaction Equation')
        #print('-----   ----------   ----------------------------------')
        #for m in range(gas.n_reactions):
        #    #print(m)
        #    print('{: 5d}   {: 10.3e}   {}'.format(m, sens1[m], gas.reaction_equation(m)))
        # use argsort to obtain an array of the *indicies* of the values of the sens1 array sorted by absolute value 
        newOrder     = np.argsort(np.abs(sens1))
        # argsort ranks from small to large, but we want large to small, so we flip this around 
        newOrder     = newOrder[::-1] 
        # make some storage variables so that we can plot the results later in a bar graph
        IndexSpecies = np.zeros(len(newOrder))
        newOrder2    = np.zeros(len(newOrder))
        sens2        = np.zeros(len(newOrder))
        reactionList = []
        # using the same method above, print the sensitivties but call the new indicies that we defined 
        print()
        print('Rxn #   k/S*dS/dk    Reaction Equation')
        print('-----   ----------   ----------------------------------')
        for ii in range(gas.n_reactions):
            print('{: 5d}   {: 10.3e}   {}'.format(newOrder[ii], sens1[newOrder[ii]], gas.reaction_equation(newOrder[ii])))
            # assign new variables values for plot use
            IndexSpecies[ii] = newOrder[ii]
            newOrder2[ii]    = ii
            sens2[ii]        = sens1[newOrder[ii]]
            reactionList.append(gas.reaction_equation(newOrder[ii]))
        # generate horizontal bar graph 
        numTopReactions = 11; # how many of the top reactions do we want to look at?
        #+++++++++++++++++++++++++++++++++++++++++++++++
        # Plotting
        #+++++++++++++++++++++++++++++++++++++++++++++++
        plt.rcdefaults()
        fig, ax = plt.subplots()
        # plot results
        ax.barh(newOrder2[0:numTopReactions-1], sens2[0:numTopReactions-1], align='center', color='black', ecolor='black')
        # make sure that every single tick on the y-axis is marked
        ax.set_yticks(np.arange(len(reactionList[0:numTopReactions-1])))
        # label these y-ticks with the reaction in question
        ax.set_yticklabels(reactionList[0:numTopReactions-1])
        # invert the y-axis so that the most sensitive reaction is at the top instead of the bottom 
        ax.invert_yaxis()
        fig.tight_layout()
        #plt.show();
        plt.savefig( str(FinalPath) + '\\SENS\\'+ str(filename) + '_PHI_' + str(phi) + '_FlameSENS.png', dpi=800);
        #plt.close()
        IndexSpeciations = [int(x) for x in IndexSpecies[0:numTopReactions-1]]
        dfXs  = pd.DataFrame({"RxnIndex":IndexSpeciations, "RxnName":reactionList[0:numTopReactions-1], "RxnSens":sens2[0:numTopReactions-1]})
        dfXs2 = pd.DataFrame({"RxnIndex":IndexSpeciations, "RxnName":reactionList[0:numTopReactions-1], "RxnSens":sorted(sens2[0:numTopReactions-1])})
        dfXs.to_csv(str(FinalPath) + '\\SENS\\'+ str(filename) + '_PHI_' + str(phi) + '_FlameSENS.dat', sep="\t", index = False)
        dfXs2.to_csv(str(FinalPath) + '\\SENS\\'+ str(filename) + '_PHI_' + str(phi) + '_FlameSENS_sorted.dat', sep="\t", index = False)
    else:
        pass
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    # Saving the DataFrame data in a *.dat format file
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    SimTime = time.time() - start
    if SimTime <= 60:
       df2 = pd.DataFrame({"File_Name":filename,"SimTime /s":round(SimTime,2),"DATE":NOW},index=[0])
    elif SimTime > 60:
       df2 = pd.DataFrame({"File_Name":filename,"SimTime /min":round((float(SimTime)/60),2),"DATE":NOW},index=[0])
    elif SimTime > 3600:
       df2 = pd.DataFrame({"File_Name":filename,"SimTime /hr":round((float(SimTime)/3600),2),"DATE":NOW},index=[0])
    df1 = pd.DataFrame({"Tinitial /K":Ti,"Pinitial /atm":Pi,"Eq.Ratio":PHI,"Exp.Result cm/s":LVB,"Sim.Result cm/s":SPEED,"% Error":ERROR,"FittedTrueSpeed cm/s and (%) Error":FITTED})
    df = df1.join(df2,lsuffix='_df1', rsuffix='_df2')
    ResultsJVB = df.to_csv(FinalPath+'\\'+mech3+'_LBV_EXP_.dat',index=False, sep="\t")
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    # Returning data
    #+++++++++++++++++++++++++++++++++++++++++++++++++
    return (ResultsJVB)
#+++++++++++++++++++++++
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
if __name__ == "__main__":
   "\t>JUST RUNNING THE CODE"
   print(half);
   mech_list = []; 
   path = "C:\\DKM\\MODELS\\" #os.getcwd()
   fname = glob.glob(path+'\\*.cti');
   # +++++++++++++++++++++++++++++++++++
   for i,mechanism in enumerate(fname):
       Input = sys.argv[1]#"InputFileNameWithFullPath"
       print("\t> Running LBV reactor for Flame speed calculation");
       print('\t> Running with mechanism: \n\t',mechanism);
       print("\t> Case: \n\t",Input)
       print(half,'\n');
       LBVCanteraSimulator(mechanism,Input);
   # +++++++++++++++++++++++++++++++++++
