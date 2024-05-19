###########################################
#  Created in 2019
#  @author: version 1.0. Sergio Martinez
###########################################
import os
import numpy as np
import pandas as pd
import os
import sys
import time
import matplotlib.pyplot as plt
import cantera as ct
import glob
import math
import datetime
import warnings
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

#from InputFileReader import InputFileReaderJSR as R
#########################################################################################
#PLOTTER FOR JSR ONLY
#++++++++++++++++++++++++++++
def plotSpecies(results,species,conditions,pressures,phi,FinalPath,Title,author,DFn,ExpTempos,filenames):
    temps=[]
    for i in results:
        if i.solution['temperature'].values[0] not in temps:
            temps.append(i.solution['temperature'].values[0])
    filename=[]
    for i in results:
        if i.mechanism not in filename:
            filename.append(i.mechanism)
    res=[]
    for i in results:
        if i.residence_time not in res:
            res.append(i.residence_time)
    pre=[]
    for i in results:
        if i.final_pressure not in pre:
            pre.append(i.final_pressure)         
            
    length=len(conditions)
    length2=len(res)
    #length3=len(pre)
    sortedresults=[results[i::length] for i in range(length)]
    for j in np.arange(len(sortedresults)):
        sortedresults[j]=[sortedresults[j][i::length2] for i in range(length2)]
    for j in np.arange(len(sortedresults)):
        for f in np.arange(len(sortedresults[j])):
            sortedresults[j][f]=[sortedresults[j][f][i::len(filename)] for i in range(len(filename))]
    #return sortedresults                           
    species2save = {}
    counter      = 0
    for w in range(2):
        plt.rcParams["font.weight"] = "bold"
        plt.tick_params(axis= 'both', direction='out', length=4, width=2, labelsize=10)
        plt.rcParams["axes.labelweight"] = "bold"
        plt.rcParams["axes.linewidth"] = "3"
        for j in np.arange(len(sortedresults)):
         for f in np.arange(len(sortedresults[j])):
             for poop in np.arange(len(species)):               
                for mech in filenames:
                  #spec=[]
                  for c in np.arange(len(filename)):                        
                     for p in np.arange(len(phi)):                        
                         specieslist=[];
                         for d in np.arange(len(temps)):                                               
                         #print(j,f,poop,c,d)
                             specieslist.append((sortedresults[j][f][c][d].solution[species[poop]].values[0]))
                         #title    = Title + " : " + author
                         #plt.title(title, fontsize=10, weight= 'bold')
                         #try:
                         Y  = pd.to_numeric(DFn[species[poop]]) 
                         #YY = (Y).sort_values(ascending=False)
                         #print(DFn.iloc[:,0])
                         #print(ExpTempos[0])
                         plt.scatter(temps,Y, marker="s", label="Exp. data: "+author, lw=4,color='black')
                         #except:
                         #    print("\t Not experimental data to plot was found, skipping...")
                         MECHA = (mech.rsplit("\\",1)[-1]).split(".cti")[0]
                         plt.plot(temps,specieslist, color="k", label=str(MECHA), lw=3)
                         plt.title(str(species[poop]).translate(subscript) + ' with residence time: ' + str(res[f])+' sec', fontsize=14, weight= 'bold')
                         #df=pd.DataFrame({'Temperature_K':temps, species[poop]+'_molefraction':specieslist})
                         species2save[str(species[poop])] = specieslist
                         #counter += 1
                         #df.to_csv(FinalPath+'\\'+str(pressures[0])+' atm'+'_species'+species[poop]+'_phi'+str(phi[p])+'_'+FileNamePlotter+'.txt',sep='\t',index=False)
                         FileNamePlotter = (((str(filename[c]).split('\\'))[-1]).split('.cti'))[0]
                         plt.xlabel('$\it{T}  /  K$', fontsize=20, weight= 'bold'); 
                         plt.ylabel('Mole fraction', fontsize=20, weight= 'bold')
                plt.legend(fontsize="medium",frameon=True)
                plt.grid(b=True, which='major', color='grey', linestyle=':'); plt.grid(b=True, which='minor', color='grey', linestyle=':')
                plt.savefig(FinalPath+'\\'+str(pressures[0])+' atm'+'_species'+species[poop]+'_phi'+str(phi[p])+'_'+FileNamePlotter+'.png',dpi=800, bbox_inches="tight")
                OutDir = (FinalPath)#+'\\'+str(pressures[0])+' atm'+'_species'+species[poop]+'_phi'+str(phi[p])+'_'+FileNamePlotter+'.png',dpi=800, bbox_inches="tight")
                #plt.savefig(FinalPath+'\\'+str(pressures[0])+' atm'+'_species'+'_phi'+str(phi[p])+'_'+FileNamePlotter+'.png',dpi=800)
                plt.close() #print(temps, specieslist)
    df     = pd.DataFrame({'Temperature_K':temps})
    Names  = list(species2save.keys())
    Values = list(species2save.values())
    #df2    = pd.DataFrame()
    for g in range(len(species2save)):
        df2 = pd.DataFrame({Names[g]:Values[g]})
        df  = pd.concat([df,df2], axis=1)
    #print(df)
    #print(df2)
    #print(df3)
    df.to_csv(FinalPath+'\\'+str(pressures[0])+' atm'+'_species'+'_phi'+str(phi[p])+'_'+FileNamePlotter+'.txt',sep='\t',index=False)    
    #df.to_csv(FinalPath+'\\'+str(pressures[0])+' atm'+'_species'+species[poop]+'_phi'+str(phi[p])+'_'+FileNamePlotter+'.txt',sep='\t',index=False)    
    #print(species2save)
    #print(len(species2save))
    return sortedresults, df, OutDir
#++++++++++++++++++++++++++++

def plotsens(sortedresults,observables,conditions,filename,temps,top10,FinalPath,Title,author):
    if not os.path.isdir(FinalPath + '\\SA\\'):
       os.makedirs(FinalPath + '\\SA\\'); 
    linetypes=['k-','b-','r-','g-','c-','m-','y-','k-.','b-.','r-.','g-.','c-.','m-.','y-.','k--','b--','r--','g--','c--','m--','y--',
               'k:','b:','r:','g:','c:','m:','y:']
    count=0
    for j in np.arange(len(sortedresults)):
        for f in np.arange(len(sortedresults[j])):
            for o in np.arange(len(observables)):
                title    = Title + " : " + author
                plt.title(title, fontsize=10, weight= 'bold')
                plt.rcParams["font.weight"] = "bold"
                plt.tick_params(axis= 'both', direction='out', length=4, width=2, labelsize=10)
                plt.rcParams["axes.labelweight"] = "bold"
                plt.rcParams["axes.linewidth"] = "3"
                fig=plt.figure()
                ax=plt.subplot(111)
                ax.set_xlabel('Temperature (K)')
                ax.set_ylabel(r'Sensitivity coefficient, $\frac{\partial\mathrm{ln}S_u^0}{\partial\mathrm{ln}k}$')
                plt.title('Top '+str(top10)+' sensitivities for observable '+observables[o])
                for c in np.arange(len(filename)):
                    
                        sens=[]
                        sens=pd.DataFrame(sortedresults[j][f][c][0].Index[1],columns=['rxn'])
                        for d in np.arange(len(temps)):
                            #print(j,f,c,d,o)
                            sens=pd.concat([sens, pd.DataFrame(sortedresults[j][f][c][d].k_sens[0],columns=[temps[d]]*len(observables)).iloc[:,o]], axis=1)
                        maxes=sens.loc[:, sens.columns != 'rxn'].max(axis=1)
                        maxes=maxes.nlargest(top10)
                        max_index=np.array(maxes.index)
                        
                        for i in max_index:
                            plt.plot(temps,sens.iloc[i][1:],linetypes[count],label=sens['rxn'][i]+', '+os.path.splitext(filename[c])[0].split('\\')[-1])
                            count=count+1
                        box = ax.get_position()
                        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                                  fancybox=False, shadow=False, ncol=1)
                plt.savefig(FinalPath+'\\'+'sensitivities_conditions'+str(j)+'_resTime'+str(f)+'_'+observables[o]+'.png',dpi=800,bbox_inches='tight')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class model_data:
    
    def __init__(self,simtype,kinetic_sens=np.array(()),physical_sens=np.array(()),Solution=pd.DataFrame(),Index=[],pIndex=[]):
        self.k_sens=kinetic_sens
        self.p_sens=physical_sens
        self.solution=Solution
        self.Index=Index
        self.pIndex=pIndex
        self.sensitivities = self.all_sensitivities()
        self.overall_index = self.sens_index()
        self.simtype=simtype
        
    # Primarily for use with JSR data, the addition operator is overloaded to combine two datasets/simulations into one array.
    # Only use this for addition of data along first axis of sensitivities    
    def __add__(self, data):
        if self.simtype==data.simtype and (self.solution.columns==data.solution.columns).all()==True:
            #If statement above checks to ensure that the two datasets come from the same type of experiment.
            if self.k_sens.any()==True and data.k_sens.any()==True:
                temp=self.k_sens.T
                temp2=data.k_sens.T
                temp=np.dstack((temp,temp2))
                self.k_sens=temp.T
        
            if self.p_sens.any()==True and data.p_sens.any()==True:
                temp=self.p_sens.T
                temp2=data.p_sens.T
                temp=np.dstack((temp,temp2))
                self.p_sens=temp.T
            
            if self.sensitivities.any()==True and data.sensitivities.any()==True:
                temp=self.sensitivities.T
                temp2=data.sensitivities.T
                temp=np.dstack((temp,temp2))
                self.sensitivities=temp.T
        
            if self.k_sens.any()==True and data.k_sens.any()==True:
                self.Index[0]=self.Index[0]+data.Index[0]
            if self.p_sens.any()==True and data.p_sens.any()==True:
                self.pIndex[0]=self.pIndex[0]+data.pIndex[0]
            if self.sensitivities.any()==True and data.sensitivities.any()==True:
                self.overall_index=self.sens_index()
            
            temp=self.solution
            temp2=data.solution
            temp = temp.append(temp2,ignore_index=True)
            
            self.solution=temp
            
            
            return self
            
        else:
            raise Exception('Simulation types must be equal in order to add model data objects together.  Simulation must also use same mechanism.')
        
        
    #The following functions take slices of sensitivity arrays eliminating one axis to create a 2d array according to observable,position,
    # or reaction.  There is a version for kinetic, physical, and all sensitivities.  Current version for position sensitivities checks
    # if there is a sensitivity at that location, and if not it will interpolate.  This portion may be re-written later.  Validity of 
    #interpolated sensitivities may be questionable    
    def species_slice_ksens(self,species_name):
        if self.k_sens.any()==True:
            position = self.Index[2].index(species_name)
            return self.k_sens[:,:,position]
        else:
            print('No kinetic sensitivities provided by simulation')
    def reaction_slice_ksens(self,reaction):
        if self.k_sens.any()==True:
            position = self.Index[1].index(reaction)
            return self.k_sens[:,position,:]
        else:
            print('No kinetic sensitivities provided by simulation')
    def x_slice_ksens(self,location):
        if self.k_sens.any()==True:
            for j in np.arange(len(self.Index[0])):
                if self.Index[0][j]==location:
                    return self.k_sens[j,:,:]
                elif j<(len(self.Index[0])-1):
                    if self.Index[0][j]<location and self.Index[0][j+1]>location:
                        result = (self.k_sens[j+1,:,:]-self.k_sens[j,:,:])*np.divide((location-self.Index[0][j]),(self.Index[0][j+1]-self.Index[0][j]))+self.k_sens[j,:,:]
                        print('Returning an interpolated result')
                        return result
                else:
                    print('Invalid grid location: check position and submit a location in range of flame solution')
        else:
            print('No kinetic sensitivities provided by simulation')
    def species_slice_psens(self,species_name):
        if self.p_sens.any()==True:
            position = self.pIndex[2].index(species_name)
            return self.p_sens[:,:,position]
        else:
            print('No physical sensitivities provided by simulation')
    def parameter_slice_psens(self,parameter):
        if self.p_sens.any()==True:
            position = self.pIndex[1].index(parameter)
            return self.p_sens[:,position,:]
        else:
            print('No physical sensitivities provided by simulation')
    def x_slice_psens(self,location):
        if self.p_sens.any()==True:
            for j in np.arange(len(self.pIndex[0])):
                if self.pIndex[0][j]==location:
                    return self.p_sens[j,:,:]
                elif j<(len(self.pIndex[0])-1):
                    if self.pIndex[0][j]<location and self.pIndex[0][j+1]>location:
                        result = (self.p_sens[j+1,:,:]-self.p_sens[j,:,:])*np.divide((location-self.pIndex[0][j]),(self.pIndex[0][j+1]-self.pIndex[0][j]))+self.p_sens[j,:,:]
                        print('Returning an interpolated result')
                        return result
                else:
                    print('Invalid grid location: check position and submit a location in range of flame solution')
        else:
            print('No physical sensitivities provided by simulation') 
    
    def species_slice_sens(self,species_name):
        if self.sensitivities.any()==True:
            position = self.overall_index[2].index(species_name)
            return self.sensitivities[:,:,position]
        else:
            print('No sensitivities provided by simulation')
    def parameter_slice_sens(self,parameter):
        if self.sensitivities.any()==True:
            position = self.overall_index[1].index(parameter)
            return self.sensitivities[:,position,:]
        else:
            print('No sensitivities provided by simulation')
    def x_slice_sens(self,location):
        if self.sensitivities.any()==True:
            for j in np.arange(len(self.overall_index[0])):
                if self.overall_index[0][j]==location:
                    return self.sensitivities[j,:,:]
                elif j<(len(self.overall_index[0])-1):
                    if self.overall_index[0][j]<location and self.overall_index[0][j+1]>location:
                        result = (self.sensitivities[j+1,:,:]-self.sensitivities[j,:,:])*np.divide((location-self.overall_index[0][j]),(self.overall_index[0][j+1]-self.overall_index[0][j]))+self.sensitivities[j,:,:]
                        print('Returning an interpolated result')
                        return result
                else:
                    print('Invalid grid location: check position and submit a location in range of flame solution')
        else:
            print('No sensitivities provided by simulation')

    def species_slice_MFExp(self,species_name):
        if self.expMolFracData.any() == True:
            if (int(np.shape(self.expMolFracData)[0])) > 1:
                position = self.MFExpIndex[1].index(species_name)
                return self.expMolFracData[:,position]
            else:
                if species_name in self.MFExpIndex[1]:
                    return self.expMolFracData
                else:
                    print('Species is not in the index')
        else:
            print('No Mole Fraction Data Provided')

            
    def sepcies_slice_MF_differences(self,species_name):
        if self.st_mf_differences.any() == True:
            if (int(np.shape(self.st_mf_differences)[0])) > 1:
                position = self.MFExpIndex[1].index(species_name)
                return self.st_mf_differences[:,position]
            else:
                if species_name in self.MFExpIndex[1]:
                    return self.st_mf_differences
                else:
                    print('That species is not in the index')
        else:
            print('No Mole Fraction Data Provided')

    def wavelength_slice_AExp(self,wavelength):
        if self.expAbsorbanceData.any() == True:
            if (int(np.shape(self.expAbsorbanceData)[0])) > 1:
                position = self.AExpIndex[1].index(wavelength)
                return self.expAbsorbanceData[:,position]
            else:
                if wavelength in self.AExpIndex[1]:
                    return self.expAbsorbanceData
                else:
                    print('That wavelength is not in the index')
        else:
            print('No Absorption Data Provided')
            
    def wavelength_slice_A_differences(self,wavelength):
        if self.st_abs_differences.any() == True :
            if (int(np.shape(self.st_abs_differences)[0])) > 1:
                position = self.AExpIndex[1].index(wavelength)
                return self.expAbsorbanceData[:,position]
            else:
                if wavelength in self.AExpIndex[1]:
                    return self.expAbsorbanceData
                else:
                    print('That wavelength is not in the index')
        else:
            print('No Absorption Data Provided')


    def Absorb_sens_slice(self,wavelength):
        if bool(self.absorbance_sens) == True:
            return self.absorbance_sens[wavelength]
        else:
            print('No Absorption Data Provided')
    def Absorb_sens_A_slice(self,wavelength):
        if bool(self.absorb_a_sens) == True:
            return self.absorb_a_sens[wavelength]
        
    def Absorb_sens_n_slice(self,wavelength):
        if bool(self.absorb_n_sens) == True:
            return self.absorb_n_sens[wavelength]
        
    def Absorb_sens_Ea_slice(self,wavelength):
        if bool(self.absorb_Ea_sens) == True:
            return self.absorb_Ea_sens[wavelength]  
          
    def reaction_slice_Absorb_sens(self,reaction,wavelength):
        if bool(self.absorbance_sens) == True:
            position = self.Index[1].index(reaction)
            return self.absorbance_sens[wavelength][:,position]
        else:
            print('No AbsoptionData Provided')
    
    def species_slice_CExp(self,species_name):
        if self.expConcentrationData.any() == True:
            if (int(np.shape(self.expConcentrationData)[0])) > 1:
                position = self.CExpInex[1].index(species_name)
                return self.expConcentrationData[:,position]
            else:
                if species_name in self.CExpInex[1]:
                    return self.expConcentrationData
                else:
                    print('That species is not in the index')
        else:
            print('No Concentration Data Provided')
    
    def species_slice_C_differences(self,species_name):
        if self.st_conc_differences.any() == True:
            if (int(np.shape(self.st_conc_differences)[0])) > 1:
                position = self.CExpInex[1].index(species_name)
                return self.st_conc_differences[:,position]
            else:
                if species_name in self.CExpInex[1]:
                    return self.st_conc_differences
                else:
                    print('That species is not in the index')
        else:
            print('No Concentration Data Provided')
    def species_slice_ksens_mappedA(self,species_name):
        if self.a_sens.any()==True:
            position = self.Index[2].index(species_name)
            return self.a_sens[:,:,position]
        else:
            print('No kinetic sensitivities provided by simulation')
            
    def reaction_slice_ksens_mappedA(self,reaction):
        if self.a_sens.any()==True:
            position = self.Index[1].index(reaction)
            return self.a_sens[:,position,:]
        else:
            print('No kinetic sensitivities provided by simulation')
    def species_slice_ksens_mappedn(self,species_name):
        if self.n_sens.any()==True:
            position = self.Index[2].index(species_name)
            return self.n_sens[:,:,position]
        else:
            print('No kinetic sensitivities provided by simulation')
            
    def reaction_slice_ksens_mappedn(self,reaction):
        if self.n_sens.any()==True:
            position = self.Index[1].index(reaction)
            return self.n_sens[:,position,:]
        else:
            print('No kinetic sensitivities provided by simulation')
    def species_slice_ksens_mappedEa(self,species_name):
        if self.e_sens.any()==True:
            position = self.Index[2].index(species_name)
            return self.e_sens[:,:,position]
        else:
            print('No kinetic sensitivities provided by simulation')
            
    def reaction_slice_ksens_mappedEa(self,reaction):
        if self.a_sens.any()==True:
            position = self.Index[1].index(reaction)
            return self.a_sens[:,position,:]
        else:
            print('No kinetic sensitivities provided by simulation')            
            
    #The following function exists to return an array of the profile of a sensitivity in one parameter with respect to an observable as
    #as a function of the independent variable (location,time, temperature)        
    def sensitivity_profile(self,observable,parameter):
        position = self.overall_index[2].index(observable)
        position_param = self.overall_index[1].index(parameter)
        sens = self.sensitivities[:,position_param,position]
        independentVar = np.asarray(self.overall_index[0])
        sens = np.flatten(sens)
        return np.vstack((independentVar,sens))
    
    #all_sensitivities returns an array which is a concatenation of physical and kinetic sensitivities available in the model
    def all_sensitivities(self):
        
        if self.p_sens.any()==False and self.k_sens.any()!=False:
            return self.k_sens
        if self.p_sens.any()!=False and self.k_sens.any()==False:
            return self.p_sens
        if self.p_sens.any()!=False and self.k_sens.any()!=False:
            return np.hstack((self.p_sens,self.k_sens))
        else:
            return np.array(())
            
    def reaction_slice_Asens(self,reaction):
        if self.a_sens.any()==True:
            position = self.Index[1].index(reaction)
            return self.a_sens[:,position,:]
        else:
            print('No A sensitivities provided by simulation')    
            
            
    def reaction_slice_nsens(self,reaction):
        if self.n_sens.any()==True:
            position = self.Index[1].index(reaction)
            return self.n_sens[:,position,:]
        else:
            print('No n sensitivities provided by simulation')
    def reaction_slice_Esens(self,reaction):
        if self.e_sens.any()==True:
            position = self.Index[1].index(reaction)
            return self.e_sens[:,position,:]
        else:
            print('No Ea sensitivities provided by simulation') 
    #sens_index returns the concatenation of the indices for the kinetic and physical sensitivities in the model.  physical
    #parameters are placed on same axis as reaction rate constants, but before the first reaction constant (eg. T,P,k1,k2,...)    
    def sens_index(self):
        if self.p_sens.any()==False and self.k_sens.any()!=False:
            return self.Index
        if self.p_sens.any()!=False and self.k_sens.any()==False:
            return self.pIndex
        if self.p_sens.any()!=False and self.k_sens.any()!=False:
            return [self.Index[0],self.pIndex[1]+self.Index[1],self.Index[2]]
        else:
            return self.Index
    def assign_flamespeed_sens(self,fsens,equations):
        self.flamespeed_sens=pd.DataFrame(data=[],index=equations)
        self.flamespeed_sens['Su']=fsens
    def plot_flamespeed_sens(self,phi,model):
        threshold = 0.02

        self.firstColumn = self.flamespeed_sens.columns[0]

        # For plotting, collect only those steps that are above the threshold
        # Otherwise, the y-axis gets crowded and illegible
        self.sensitivitiesSubset = self.flamespeed_sens[self.flamespeed_sens[self.firstColumn].abs() > threshold]
        #print(sensitivitiesSubset)
        self.indicesMeetingThreshold = self.sensitivitiesSubset[self.firstColumn].abs().sort_values(ascending=False).index
        self.sensitivitiesSubset.loc[self.indicesMeetingThreshold].plot.barh(title="Sensitivities for "+model+" at phi="+str(phi),
                                                          legend=None)
        plt.gca().invert_yaxis()

        plt.rcParams.update({'axes.labelsize': 20})
        plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');

        # Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
        plt.savefig('sensitivityPlot', dpi=600)
    def add_mechanism(self,mech):
        self.mechanism=mech
    def add_yaml(self,yaml):
        self.yaml_input=yaml
    def assign_phi(self,phi):
        self.phi=phi
    def add_fuel(self,fuel):
        self.fuel=fuel
    def add_exp_shocktube_MolFracs(self,molFracData= np.array(()),differences= np.array(()), index=[]):
        self.expMolFracData=molFracData
        self.st_mf_differences=differences
        self.MFExpIndex= index
    def add_exp_shocktube_absorbance(self,absorbanceData = np.array(()),differences = np.array(()),absorbanceProfiles = {},absorbanceSens = {}, index = [],absorbancePSens =[]):
        self.expAbsorbanceData=absorbanceData
        self.st_abs_differences=differences
        self.absorbance_profiles=absorbanceProfiles
        self.absorbance_sens = absorbanceSens
        self.AExpIndex = index
        self.absorbance_p_sens = absorbancePSens
    def add_exp_shocktube_concentrations(self,concData = np.array(()),differences = np.array(()),index = []):
        self.expConcentrationData=concData
        self.st_conc_differences=differences
        self.CExpInex = index
    def add_S_Matrix(self, SMatrix = np.array(())):
        self.S_matrix = SMatrix
    def add_Y_Matrix(self,YMatrix = np.array(())):
        self.Y_matrix = YMatrix
    def assign_ksens_mappings(self, asens = np.array(()),nsens= np.array(()),easens= np.array(())):
        self.a_sens=asens
        self.n_sens=nsens
        self.e_sens=easens
    def absorbance_sens_mappings(self,asens = {}, nsens = {}, easens = {}):
        self.absorb_a_sens = asens
        self.absorb_n_sens = nsens
        self.absorb_Ea_sens = easens
    def add_parameter_arrays(self,molecularParameterArray = [],molecularParameterArrayCombinedReaction = [], fullArray = np.array(())):
        self.molecular_parameter_sens = molecularParameterArray
        self.molecular_parameter_sens_combined_reaction = molecularParameterArrayCombinedReaction
        self.full_array_parameter_sens = fullArray
    def write_model_to_file(self,filename):
        items=self.__dict__.keys()
        #print(dir(self))
        with open(filename,'w') as f:
            for j in items:
                f.write(j+'\n')
                if j=='simtype':
                    f.write(getattr(self,j)+'\n')
                elif j=='Index' and self.Index!=[]:
                    for k in np.arange(len(self.Index)):
                        if k==0:
                            f.write('    z_T_t\n')
                        elif k==1:
                            f.write('    reactions\n')
                        elif k==2:
                            f.write('    observables\n')
                        for val in self.Index[k]:
                            f.write(str(val)+'\n')
                elif j=='overall_index' and self.overall_index!=[]:
                    for k in np.arange(len(self.overall_index)):
                        if k==0:
                            f.write('    z_T_t\n')
                        elif k==1:
                            f.write('    reactions and physical_params\n')
                        elif k==2:
                            f.write('    observables\n')
                        for val in self.Index[k]:
                            f.write(str(val)+'\n')
                elif j=='final_pressure':
                    f.write(str(getattr(self,j))+'\n')
                elif j=='pIndex' and self.pIndex!=[]:
                    for k in np.arange(len(self.pIndex)):
                        if k==0:
                            f.write('    z_T_t\n')
                        elif k==1:
                            f.write('    physical_params\n')
                        elif k==2:
                            f.write('    observables\n')
                        for val in self.Index[k]:
                            f.write(str(val)+'\n')
                elif j=='mechanism':
                    f.write(self.mechanism+'\n')
                elif j=='yaml_input':
                    f.write(self.yaml_input+'\n')
                elif j=='phi':
                    f.write(str(self.phi)+'\n')
                elif j=='fuel':
                    f.write(self.fuel+'\n')
                elif j=='solver_time':
                    f.write(str(self.solver_time)+'\n')
                #f.writelines(getattr(self,j)+['\n'])
    def add_interpolated_temps_for_physical_sens(self,listOfTemperatureArrays = []):
        self.temps_for_physical_sens = listOfTemperatureArrays
    def add_interpolated_pressure_for_physical_sens(self,listOfPressureArrays = []):
        self.pressure_for_physical_sens = listOfPressureArrays
    def add_interpolated_concentration_for_physical_sens(self,listOfConcentrationDataFrames = []):
        self.concentration_for_physical_sens = listOfConcentrationDataFrames
    def add_psens_array(self,fullArray = np.array(())):
        self.full_pSens_array = fullArray
    def add_target_values(self,initalTemperature,initialPressure,initialMoleFractions,initialKvalues=[],initialMolecularProperties=[]):
        self.initial_temp = initalTemperature
        self.initial_pressure = initialPressure
        self.initial_MF = initialMoleFractions
        self.target_rate_constants = initialKvalues
        self.target_molecular_props = initialMolecularProperties
        
    def add_final_pressure(self,pressure):
        self.final_pressure=pressure
    def add_solver_time(self,solvertime):
        self.solver_time=solvertime
    def add_residence_time(self,res):
        self.residence_time=res
    def add_temperature_sensitivity(self,sens):
        self.temp_sens=np.array(sens)
    def add_forward_rates(self,rates):
        self.forward_rates=rates
    def add_reverse_rates(self,rates):
        self.reverse_rates=rates
    def add_net_rates_of_progress(self,rates):
        self.net_rates_of_progress=rates
    def tags(self,tags):
        self.tags=tags
#This function is for pre-mixed burner flame
# *******************************************************************************

def InputFileReaderJSR(FileNamePath,OutputPath):
    FILENAME = FileNamePath.split('\\')
    INPUTFILENAME = FILENAME[-1]; 
    FileJSR= pd.read_csv(FileNamePath, names=['col'],sep='%s',engine='python');
    FolderName = FILENAME[3:-1]; 
    InputFileNames = (INPUTFILENAME.split(".JSR_EXP"))[0]; 
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
    try:
        REACS = (FileJSR[FileJSR.col.str.contains("REAC:")]).col.str.split(':')
    except:
        FileJSR = pd.read_fwf(FileNamePath,names=['col']);
        REACS  = (FileJSR[FileJSR.col.str.contains("REAC:")]).col.str.split(':')
    REACTANTS = []
    for i in range(len(REACS)):
        REAC1 = REACS.iloc[i][1]
        REAC1 = ":".join(REAC1.split())
        REACTANTS.append(REAC1)
    REACTANTS = str(",".join(REACTANTS))
    AUTHOR = ((FileJSR[FileJSR.col.str.contains("AUTHOR:")]
               ).iloc[0].str.split(':'))[0][1]
    JOURNAL = str(FileJSR.iloc[1,:])#((FileJSR[FileJSR.col.str.contains("JOURNAL:")]
#                ).iloc[0].str.split(':'))[0][1]
    PDF = ((FileJSR[FileJSR.col.str.contains("PDF:")]
            ).iloc[0].str.split(':'))[0][1]
    REACTOR = ((FileJSR[FileJSR.col.str.contains("REACTOR:")]
                ).iloc[0].str.split(':'))[0][1]
    DATA = ((FileJSR[FileJSR.col.str.contains("DATA:")]
             ).iloc[0].str.split(':'))[0][1]
    PRESSURE = ((FileJSR[FileJSR.col.str.contains("PRESSURE:")]
                 ).iloc[0].str.split(':'))[0][1]
    PRESSURE = str(PRESSURE)
    PHI = ((FileJSR[FileJSR.col.str.contains("PHI:")]
            ).iloc[0].str.split(':'))[0][1]
    E_EQUATION = ((FileJSR[FileJSR.col.str.contains("ENERGY_EQUATION:")]
            ).iloc[0].str.split(':'))[0][1]
    try:
        T_INITIAL = ((FileJSR[FileJSR.col.str.contains("T_INITIAL:")]
                      ).iloc[0].str.split(':'))[0][1]
        P_INITIAL = ((FileJSR[FileJSR.col.str.contains("P_INITIAL:")]
                        ).iloc[0].str.split(':'))[0][1]
        ReactorsDict["T_INITIAL"] = T_INITIAL
        ReactorsDict["P_INITIAL"] = P_INITIAL
    except:
        ReactorsDict["T_INITIAL"] = 0
        ReactorsDict["P_INITIAL"] = 0
    TEMPERATURE = ((FileJSR[FileJSR.col.str.contains(
        "TEMPERATURE:")]).iloc[0].str.split(':'))[0][1]
    TIME = ((FileJSR[FileJSR.col.str.contains("TIME:")]
             ).iloc[0].str.split(':'))[0][1]
    DSpeciesJSR     = FileJSR[FileJSR.col.str.contains("!") == False]
    DSJSR           = DSpeciesJSR[DSpeciesJSR.col.str.contains("TAU") == True]
    DataEspeciesJSR = pd.DataFrame((DSJSR.col.str.split()).tolist())
    DropingStuff    = FileJSR[FileJSR.col.str.contains("!") == False]
    DropingStuff2   = DropingStuff[DropingStuff.col.str.contains("TAU") == False]
    dataJSR         = pd.DataFrame((DropingStuff2.col.str.split()).tolist())
    NamesDict["AUTHOR"]             = AUTHOR
    NamesDict["JOURNAL"]            = JOURNAL
    NamesDict["PDF"]                = PDF
    NamesDict["InpFileName"]        = INPUTFILENAME
    NamesDict["FolderDir"]          = FolderDir
    ReactorsDict["REACTOR"]         = REACTOR
    ReactorsDict["DATA"]            = DATA
    ReactorsDict["REACTANTS"]       = REACTANTS
    ReactorsDict["PHI"]             = PHI
    ReactorsDict["ENERGY_EQUATION"] = E_EQUATION
    UnitsDict["PRESSURE"]           = PRESSURE
    UnitsDict["TEMPERATURE"]        = TEMPERATURE
    UnitsDict["TIME"]               = TIME
    return (NamesDict, ReactorsDict, UnitsDict,DataEspeciesJSR,dataJSR,FinalPath)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#This function is for pre-mixed burner flame
def burner_flame(gas,grid,mdot,data=pd.DataFrame(columns=['z','T']),kinetic_sens=0,physical_sens=0,observables=[],physical_params=['T','P'],energycon=False,soret=True):
    #when energycon is off treat flame as burner stabilized, with known T-profile
    simtype = 'burner flame'
    baseConditions=gas.TPX
    tol_ss = [1.0e-5, 1.0e-13]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-4, 1.0e-10]  # [rtol atol] for time stepping
    loglevel = 1  # amount of diagnostic output (0 to 5)
   
    f = ct.BurnerFlame(gas, width=grid)
    f.burner.mdot = mdot
    
    f.set_initial_guess()
    # read temperature vs. position data from a file.
    # The file is assumed to have one z, T pair per line, separated by a comma.
    if data.empty==False and energycon==False:
        zloc=data['z']
        tvalues=data['T']  
        zloc /= max(zloc)

        #print(tvalues)

        f.flame.set_fixed_temp_profile(zloc, tvalues)  #sets a fixed temperature profile for the flame simulation.  May come from a measurement. Requires no energy conservation
    elif data.empty==False and energycon!='off':
        raise Exception('User has supplied fixed temperature dataset but energy conservation is not off.  Remove dataset or turn energy conservation off')
        
        
    f.flame.set_steady_tolerances(default=tol_ss)  #Set steady tolerances
    f.flame.set_transient_tolerances(default=tol_ts) #Set transient tolerances
    f.show_solution()

    f.energy_enabled = energycon  #This must be set to false for a burner stabilized flame with known T-profile

    
    
    f.transport_model = 'Multi'   #Sets to multicomponent transport for simulation.  Needs to be set this way to use Soret effect
    f.set_max_jac_age(10, 10)       #Age limits on Jacobian-leave as is for best results
    f.solve(loglevel, refine_grid=False)  #Solve for initial estimate without grid refinement
    
    f.soret_enabled = soret          #Enable Soret effect.  Remember transport must be set to multi.  Mix causes failure

    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.5)  #Establishes refinement criteria for grid

    #print('mixture-averaged flamespeed = ', f.u[0])

    f.transport_model = 'Multi'         #This block solves problem again with grid refinement on
    f.solve(loglevel, refine_grid=True)
    f.show_solution()
    print('multicomponent flamespeed = ', f.u[0])
    

    
    #solution = f
    ##Begin section to calculate sensitivities
    dk = 0.010
    solution = f.X
    if kinetic_sens==1 and bool(observables):
        
        #Calculate kinetic sensitivities
        sensIndex = [f.grid.tolist(),gas.reaction_equations(),observables]
        
        S = np.zeros((len(f.grid),gas.n_reactions,len(observables)))
        #print(solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1])
        #a=solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1]
        for m in range(gas.n_reactions):
            gas.set_multiplier(1.0)
            gas.set_multiplier(1+dk,m)
            f.solve(loglevel=1,refine_grid=False)
            for i in np.arange(len(observables)):
                for k in np.arange(len(f.grid)):                    
                    S[k,m,i]=f.X[f.flame.component_index(observables[i])-4,k]-solution[f.flame.component_index(observables[i])-4,k]
                    #print(solution.X[solution.flame.component_index(observables[i])-4,k])
                    #print(f.X[f.flame.component_index(observables[i])-4,k])
                    S[k,m,i]=np.divide(S[k,m,i],solution[f.flame.component_index(observables[i])-4,k])             
                    S[k,m,i]=np.divide(S[k,m,i],dk)
                    
                    
    
    if physical_sens==1 and bool(observables):
        #Calculate physical sensitivities
        gas.set_multiplier(1.0)
        
        psensIndex = [f.grid.tolist(),physical_params,observables]
        pS = np.zeros((len(f.grid),len(physical_params),len(observables)))
        for m in range(len(physical_params)):
            gas.TPX=baseConditions
            if physical_params[m]=='T':
                gas.TPX=baseConditions[0]+dk,baseConditions[1],baseConditions[2]
            elif physical_params[m]=='P':
                gas.TPX=baseConditions[0],baseConditions[1]+dk,baseConditions[2]
            f.solve(loglevel=1,refine_grid=False)
            for i in np.arange(len(observables)):
                for k in np.arange(len(f.grid)):
                    pS[k,m,i] =np.log10(solution[f.flame.component_index(observables[i])-4,k])-np.log10(f.X[f.flame.component_index(observables[i])-4,k])
                    pS[k,m,i] = np.divide(pS[k,m,i],np.log10(dk))
                    
                    
                    
    
    elif kinetic_sens==1 and bool(observables)==False:
        raise Exception('Please supply a list of observables in order to run kinetic sensitivity analysis')
    elif physical_sens==1 and bool(observables)==False:
        raise Exception('Please supply a list of observables in order to run physical sensitivity analysis')
    gas.set_multiplier(1.0)
    gas.TP = baseConditions[0],baseConditions[1]  
    f.solve(loglevel=1,refine_grid=False)
    solution=pd.DataFrame(columns=f.flame.component_names)
    #for i in f.flame.component_names:
        #solution[i]=f.solution(i)
    for i in np.arange(len(f.flame.component_names)):
        if i<=3:
            solution[f.flame.component_names[i]]=f.solution(i)
        else:
            solution[f.flame.component_names[i]]=f.X[i-4,:]
       
    if kinetic_sens==1 and bool(observables) and physical_sens!=1:
        results = model_data(simtype,kinetic_sens=S,Solution=solution,Index=sensIndex)        
        return results
    elif physical_sens==1 and bool(observables) and kinetic_sens!=1:
        results = model_data(simtype,Solution=solution,pIndex=psensIndex,physical_sens=pS)
        return results
    elif kinetic_sens==1 and physical_sens==1 and bool(observables):
        results = model_data(simtype,Solution=solution,pIndex=psensIndex,Index=sensIndex,physical_sens=pS,kinetic_sens=S)
        return results
    elif kinetic_sens!=1 and physical_sens!=1:
        Index = [f.grid.tolist()] 
        results = model_data(simtype,Solution=solution,Index=Index)
        return results
    else:
        print('Something went wrong with the parameters given.  Kinetic and physical sens may be set to either 0 or 1, and require a list of observables')

def JSR_steadystate(gas,resTime,volume,kinetic_sens=0,physical_sens=0,observables=[],physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=1000):
    # Inlet gas conditions are passed into function in the "gas" object, which is a cantera object
    # Reactor parameters passed into function as resTime and volume.  Residence time and volume of JSR
    #kinetic sens and physical sens are optional parameters which are set to zero by default.  Set them to 1 to
    #calculate sensitivity based on kinetic or physical parameters.  If these are set to 1 you must pass
    #an array of all observables to calculate sensitivities for
    
    simtype='jsr'
    reactorPressure=gas.P
    # This is the "conductance" of the pressure valve and will determine its efficiency in 
    # holding the reactor pressure to the desired conditions. It is an optional parameter
    pressureValveCoefficient=pressureValveCoefficient

    # This parameter will allow you to decide if the valve's conductance is acceptable. If there
    # is a pressure rise in the reactor beyond this tolerance, you will get a warning
    maxPressureRiseAllowed = 0.001

    # Simulation termination criterion
    #maxSimulationTime = maxsimulationTime # seconds.  An optional parameter
    fuelAirMixtureTank = ct.Reservoir(gas)
    exhaust = ct.Reservoir(gas)

    stirredReactor = ct.IdealGasReactor(gas, energy=energycon, volume=volume)
    
    
        
    massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,downstream=stirredReactor,mdot=stirredReactor.mass/resTime)

    pressureRegulator = ct.Valve(upstream=stirredReactor,downstream=exhaust,K=pressureValveCoefficient)

    reactorNetwork = ct.ReactorNet([stirredReactor])
    
    #This block adds kinetic sensitivity parameters for all reactions if desired.  
    if kinetic_sens==1 and bool(observables):
        for i in range(gas.n_reactions):
            stirredReactor.add_sensitivity_reaction(i)
            
    if kinetic_sens==1 and bool(observables)==False:
        print('Please supply a non-empty list of observables for sensitivity analysis or set kinetic_sens=0')
        
    #if physical_sens==1 and bool(observables):
        #print('Placeholder')
    if physical_sens==1 and bool(observables)==False:
        print('Please supply a non-empty list of observables for sensitivity analysis or set physical_sens=0')

    # now compile a list of all variables for which we will store data
    columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    columnNames = ['pressure'] + columnNames

    # use the above list to create a DataFrame
    timeHistory = pd.DataFrame(columns=columnNames)

    # Start the stopwatch
    tic = time.time()
    
    
    #Names=[]
    #for l in np.arange(stirredReactor.n_vars):
        #Names.append(stirredReactor.component_name(l))
    #global b
    
    #Names = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    #state = np.hstack([stirredReactor.mass, 
                   #stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
    #print(state)
    #b=pd.DataFrame(data=state).transpose()
    #b.columns=Names
    
    #Establish a matrix to hold sensitivities for kinetic parameters, along with tolerances
    if kinetic_sens==1 and bool(observables):
        #senscolumnNames = ['Reaction']+observables     
        senscolumnNames = observables
        #sensArray = pd.DataFrame(columns=senscolumnNames)
        senstempArray = np.zeros((gas.n_reactions,len(observables)))
        reactorNetwork.rtol_sensitivity = 1.0e-6
        reactorNetwork.atol_sensitivity = 1.0e-6
        
    dk=0.01    
        
    if physical_sens==1 and bool(observables):
        #psenscolumnNames = ['Parameter'] + observables
        #psensArray = pd.DataFrame(columns=senscolumnNames)
        pIndex=[[gas.T],physical_params,observables]
        psenstempArray = np.zeros((len(observables),len(physical_params)))
        tempSol=[]
        conditions=gas.TPX
        for i in np.arange(len(physical_params)):
            if physical_params[i]=='T':
                gas.TPX=conditions[0]+dk,conditions[1],conditions[2]
            if physical_params[i]=='P':
                gas.TPX=conditions[0],conditions[1]+dk,conditions[2]
            
            tempMixTank=ct.Reservoir(gas)
            tempExhaust = ct.Reservoir(gas)
            tempReactor=ct.IdealGasReactor(gas,energy=energycon,volume=volume)
            tempMassFlowCt=ct.MassFlowController(upstream=tempMixTank,downstream=tempReactor,mdot=tempReactor.mass/resTime)
            tempPresReg=ct.Valve(upstream=tempReactor,downstream=tempExhaust,K=pressureValveCoefficient)
            tempNetwork=ct.ReactorNet([tempReactor])
            tempNetwork.advance_to_steady_state()
            tempSol.append(tempReactor.get_state())
            gas.TPX=conditions
        
    reactorNetwork.advance_to_steady_state()
    final_pressure=stirredReactor.thermo.P
    global b
    b=reactorNetwork.sensitivities()
    if kinetic_sens==1 and bool(observables):
        for k in np.arange(len(observables)):
            for j in np.arange(gas.n_reactions):
                try:
                    senstempArray[j,k]=reactorNetwork.sensitivity(observables[k],j)
                except:
                    senstempArray[j,k]=-1
                
        #sensArray['Reaction']=gas.reaction_equations()        
        #sensArray[observables]=senstempArray.T   
        #temp = sensArray.as_matrix()
        kIndex = [[gas.T],gas.reaction_equations(),senscolumnNames]
        
        
        
    if physical_sens==1 and bool(observables):
        for k in np.arange(len(observables)):
            for j in np.arange(len(['T','P'])):
                psenstempArray[j,k]=np.log10(stirredReactor.get_state()[stirredReactor.component_index(observables[k])])-np.log10(tempSol[j][stirredReactor.component_index(observables[k])])
                psenstempArray[j,k]=np.divide(psenstempArray[j,k],np.log10(dk))
    
    #state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                   #stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X]) 
        
        # Stop the stopwatch
    toc = time.time()

    print('Simulation Took {:3.2f}s to compute'.format(toc-tic))
    columnNames = []
    #Store solution to a solution array
    #for l in np.arange(stirredReactor.n_vars):
        #columnNames.append(stirredReactor.component_name(l))
    columnNames=[stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    #state=stirredReactor.get_state()
    state=np.hstack([stirredReactor.mass, 
                   stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
    data=pd.DataFrame(state).transpose()
    data.columns=columnNames
       
    # We now check to see if the pressure rise during the simulation, a.k.a the pressure valve
    # was okay
    pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
    if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
        print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
        
    if kinetic_sens==1 and bool(observables) and physical_sens!=1:
        modelData=model_data(simtype,kinetic_sens=np.expand_dims(senstempArray,axis=0),Solution=data,Index=kIndex)
        modelData.add_final_pressure(final_pressure)
        return modelData
    if physical_sens==1 and bool(observables) and kinetic_sens!=1:
        
        modelData=model_data(simtype,physical_sens=np.expand_dims(psenstempArray,axis=0),Solution=data,pIndex=pIndex)
        modelData.add_final_pressure(final_pressure)
        return modelData
    if physical_sens==1 and bool(observables) and kinetic_sens==1:
        
        modelData=model_data(simtype,physical_sens=np.expand_dims(psenstempArray,axis=0),kinetic_sens=np.expand_dims(senstempArray,axis=0),Solution=data,Index=kIndex,pIndex=pIndex)
        modelData.add_final_pressure(final_pressure)
        return modelData
    else:
        modelData=model_data(simtype,Solution=data)
        modelData.add_final_pressure(final_pressure)
        return modelData

def multiTemp(cti,gas,Temps,f,kinetic_sens=0,physical_sens=0,observables=[],physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=1000):
    gas.TPX=Temps[0],f['pressure']*ct.one_atm,f['conditions']
    
    
    model=JSR_steadystate(gas,f['residenceTime'],f['reactorVolume'],kinetic_sens=kinetic_sens,physical_sens=physical_sens,observables=observables,physical_params=physical_params,energycon=energycon,pressureValveCoefficient=pressureValveCoefficient) 
    for i in np.arange(len(Temps)):
        if i!=0:
            try:
                
                
                gas.TPX = Temps[i],f['pressure']*ct.one_atm,f['conditions']
                solutionObject = JSR_steadystate(gas,f['residenceTime'],f['reactorVolume'],kinetic_sens=kinetic_sens,physical_sens=physical_sens,observables=observables,physical_params=physical_params,energycon=energycon,pressureValveCoefficient=pressureValveCoefficient)
            
                model=solutionObject+model
            except:
                print('Simulation at '+str(Temps[i])+' K failed.')
                pass
    return model

def JSR_steadystate2(gas,resTime,volume,kinetic_sens=0,physical_sens=0,observables=[],physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=1000,initial_conditions_gas=0,tempsens=0):
    # Inlet gas conditions are passed into function in the "gas" object, which is a cantera object
    # Reactor parameters passed into function as resTime and volume.  Residence time and volume of JSR
    #kinetic sens and physical sens are optional parameters which are set to zero by default.  Set them to 1 to
    #calculate sensitivity based on kinetic or physical parameters.  If these are set to 1 you must pass
    #an array of all observables to calculate sensitivities for
    
    simtype='jsr'
    reactorPressure=gas.P
    # This is the "conductance" of the pressure valve and will determine its efficiency in 
    # holding the reactor pressure to the desired conditions. It is an optional parameter
    pressureValveCoefficient=pressureValveCoefficient

    # This parameter will allow you to decide if the valve's conductance is acceptable. If there
    # is a pressure rise in the reactor beyond this tolerance, you will get a warning
    maxPressureRiseAllowed = 0.001

    # Simulation termination criterion
    #maxSimulationTime = maxsimulationTime # seconds.  An optional parameter
    fuelAirMixtureTank = ct.Reservoir(gas)
    exhaust = ct.Reservoir(gas)
    if initial_conditions_gas==0:
        stirredReactor = ct.IdealGasReactor(gas, energy=energycon, volume=volume)
        mdot=stirredReactor.mass/resTime
    else:
        stirredReactor = ct.IdealGasReactor(initial_conditions_gas,energy=energycon,volume=volume)
        dummyReactor = ct.IdealGasReactor(gas,energy=energycon,volume=volume)
        mdot=dummyReactor.mass/resTime
    
    
        
    massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,downstream=stirredReactor,mdot=mdot)

    pressureRegulator = ct.Valve(upstream=stirredReactor,downstream=exhaust,K=pressureValveCoefficient)

    reactorNetwork = ct.ReactorNet([stirredReactor])

    #This block adds kinetic sensitivity parameters for all reactions if desired.  
    if kinetic_sens==1 and bool(observables):
        for i in range(gas.n_reactions):
            stirredReactor.add_sensitivity_reaction(i)
            
    if kinetic_sens==1 and bool(observables)==False:
        print('Please supply a non-empty list of observables for sensitivity analysis or set kinetic_sens=0')
        
    #if physical_sens==1 and bool(observables):
        #print('Placeholder')
    if physical_sens==1 and bool(observables)==False:
        print('Please supply a non-empty list of observables for sensitivity analysis or set physical_sens=0')

    # now compile a list of all variables for which we will store data
    columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    columnNames = ['pressure'] + columnNames

    # use the above list to create a DataFrame
    timeHistory = pd.DataFrame(columns=columnNames)

    # Start the stopwatch
    tic = time.time()


    #print('Simulation Took {:3.2f}s to compute'.format(toc-tic))

    #Names=[]
    #for l in np.arange(stirredReactor.n_vars):
        #Names.append(stirredReactor.component_name(l))
    #global b
    
    #Names = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    #state = np.hstack([stirredReactor.mass, 
                   #stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
    #print(state)
    #b=pd.DataFrame(data=state).transpose()
    #b.columns=Names
    
    #Establish a matrix to hold sensitivities for kinetic parameters, along with tolerances
    if kinetic_sens==1 and bool(observables):
        #senscolumnNames = ['Reaction']+observables     
        senscolumnNames = observables
        #sensArray = pd.DataFrame(columns=senscolumnNames)
        senstempArray = np.zeros((gas.n_reactions,len(observables)))
        reactorNetwork.rtol_sensitivity = 1.0e-12
        reactorNetwork.atol_sensitivity = 1.0e-12
        
    dk=0.01    
        
    if physical_sens==1 and bool(observables):

        #psenscolumnNames = ['Parameter'] + observables
        #psensArray = pd.DataFrame(columns=senscolumnNames)
        pIndex=[[gas.T],physical_params,observables]
        psenstempArray = np.zeros((len(observables),len(physical_params)))
        tempSol=[]
        conditions=gas.TPX
        for i in np.arange(len(physical_params)):
            if physical_params[i]=='T':
                gas.TPX=conditions[0]+dk,conditions[1],conditions[2]
            if physical_params[i]=='P':
                gas.TPX=conditions[0],conditions[1]+dk,conditions[2]
            
            tempMixTank=ct.Reservoir(gas)
            tempExhaust = ct.Reservoir(gas)
            tempReactor=ct.IdealGasReactor(gas,energy=energycon,volume=volume)
            tempMassFlowCt=ct.MassFlowController(upstream=tempMixTank,downstream=tempReactor,mdot=tempReactor.mass/resTime)
            tempPresReg=ct.Valve(upstream=tempReactor,downstream=tempExhaust,K=pressureValveCoefficient)
            tempNetwork=ct.ReactorNet([tempReactor])
            tempNetwork.advance_to_steady_state()
            tempSol.append(tempReactor.get_state())
            gas.TPX=conditions
    #reactorNetwork.rtol=1e-9
    #resid=reactorNetwork.advance_to_steady_state(return_residuals=True)
    #print(resid)
    t = 0
    while t < maxsimulationTime:
        t = reactorNetwork.step()
    toc = time.time()
    print("\t> Reactor Network rtol:\n\t> {0}".format(reactorNetwork.rtol))
    final_pressure=stirredReactor.thermo.P
    global b
    b=reactorNetwork.sensitivities()
    if kinetic_sens==1 and bool(observables):
        for k in np.arange(len(observables)):
            for j in np.arange(gas.n_reactions):
                try:
                    senstempArray[j,k]=reactorNetwork.sensitivity(observables[k],j)
                except:
                    senstempArray[j,k]=-1
                
        #sensArray['Reaction']=gas.reaction_equations()        
        #sensArray[observables]=senstempArray.T   
        #temp = sensArray.as_matrix()
        kIndex = [[gas.T],gas.reaction_equations(),senscolumnNames]
        
        
        
    if physical_sens==1 and bool(observables):
        for k in np.arange(len(observables)):
            for j in np.arange(len(['T','P'])):
                psenstempArray[j,k]=np.log10(stirredReactor.get_state()[stirredReactor.component_index(observables[k])])-np.log10(tempSol[j][stirredReactor.component_index(observables[k])])
                psenstempArray[j,k]=np.divide(psenstempArray[j,k],np.log10(dk))
    
    #state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                   #stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X]) 
        
        # Stop the stopwatch
    toc = time.time()
    print('\t> Simulation at T={}K took {:3.2f}s to compute'.format(gas.T, toc-tic))

    columnNames = []
    #Store solution to a solution array
    #for l in np.arange(stirredReactor.n_vars):
        #columnNames.append(stirredReactor.component_name(l))
    columnNames=[stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    #state=stirredReactor.get_state()
    state=np.hstack([stirredReactor.mass, 
                   stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
    data=pd.DataFrame(state).transpose()
    data.columns=columnNames
       
    # We now check to see if the pressure rise during the simulation, a.k.a the pressure valve
    # was okay
    pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
    if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
        print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
        
    if kinetic_sens==1 and bool(observables) and physical_sens!=1:
        modelData=model_data(simtype,kinetic_sens=np.expand_dims(senstempArray,axis=0),Solution=data,Index=kIndex)
        modelData.add_final_pressure(final_pressure)
        modelData.add_solver_time(toc-tic)
        return modelData
    if physical_sens==1 and bool(observables) and kinetic_sens!=1:
        
        modelData=model_data(simtype,physical_sens=np.expand_dims(psenstempArray,axis=0),Solution=data,pIndex=pIndex)
        modelData.add_final_pressure(final_pressure)
        modelData.add_solver_time(toc-tic)
        return modelData
    if physical_sens==1 and bool(observables) and kinetic_sens==1:
        
        modelData=model_data(simtype,physical_sens=np.expand_dims(psenstempArray,axis=0),kinetic_sens=np.expand_dims(senstempArray,axis=0),Solution=data,Index=kIndex,pIndex=pIndex)
        modelData.add_final_pressure(final_pressure)
        modelData.add_solver_time(toc-tic)
        return modelData
    else:
        modelData=model_data(simtype,Solution=data)
        modelData.add_final_pressure(final_pressure)
        modelData.add_solver_time(toc-tic)
        return modelData

def multiTemp2(cti,gas,Temps,f,kinetic_sens=0,physical_sens=0,observables=[],physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=1000,initial_condition_gas=0,tempsens=0):
    gas.TPX=Temps[0],f['pressure']*ct.one_atm,f['conditions']
    
    
    model=JSR_steadystate2(gas,f['residenceTime'],f['reactorVolume'],kinetic_sens=kinetic_sens,physical_sens=physical_sens,observables=observables,physical_params=physical_params,energycon=energycon,pressureValveCoefficient=pressureValveCoefficient,initial_conditions_gas=initial_condition_gas,tempsens=tempsens) 
    for i in np.arange(len(Temps)):
        if i!=0:
            try:
                gas.TPX = Temps[i],f['pressure']*ct.one_atm,f['conditions']
                solutionObject = JSR_steadystate2(gas,f['residenceTime'],f['reactorVolume'],kinetic_sens=kinetic_sens,physical_sens=physical_sens,observables=observables,physical_params=physical_params,energycon=energycon,pressureValveCoefficient=pressureValveCoefficient,initial_condition_gass=initial_condition_gas,tempsens=tempsens)
            
                model=solutionObject+model
            except:
                print('Simulation at '+str(Temps[i])+' K failed.')
                pass
    return model

#+++++++++++++++++++++++++++++++++++++
# end of the definition of model_data
#+++++++++++++++++++++++++++++++++++++
#########################################################################################
# definition for SOLVER class::::
#++++++++++++++++++++++++++++++++

def JSR(filenames,Temps,pressures,resTime,conditions,sens,volume,observables):
    expConditions=[]
    for T in Temps:
        for filename in filenames:
            for P in pressures:
                for r in resTime:
                    for X in conditions:
                        expConditions.append([T,P,filename,r,X])
    tempresults=[]
    results=[]
    for i in expConditions:
        gas=ct.Solution(i[2])
        f={'residenceTime':i[3],'reactorVolume':volume,'pressure':i[1],'conditions':i[4]}
        results.append(multiTemp2(i[2],gas,[i[0]],f,kinetic_sens=0,physical_sens=0,observables=observables,physical_params=['T','P'],
        energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=10000))
        results[-1].add_mechanism(i[2])
        #results[-1].add_fuel(fuel)
        #results[-1].assign_phi(phi)
        results[-1].add_residence_time(i[3])
        print("\t> Final pressure:\n\t> {0} atm".format(round((results[-1].final_pressure)/101325,2)))
    return (results)
#++++++++++++++++++++++++++++++++++++
# End of the definition for SOLVER
#++++++++++++++++++++++++++++++++++++

#Change inputs here
def JSRCanteraSimulator(NamesDict, ReactorsDict, UnitsDict, DataEspeciesJSR, TP, FinalPath, MechFile, sa):

    ct.suppress_thermo_warnings()
    Mech = [MechFile]
    start = time.time();
    #print("NamesDict: ", NamesDict)
    #print("ReactorsDict: ", ReactorsDict)
    #print("UnitsDict: ", UnitsDict)
    #print("DataEspeciesJSR: ", DataEspeciesJSR)
    #print("TP: ", TP)
    #print("FinalPath: ", FinalPath)
    #print("MechFile: ", MechFile)
    #print("sa: ", sa)
    TINI = float(ReactorsDict["T_INITIAL"]);
    PINI = float(ReactorsDict["P_INITIAL"]);
    #print(TINI,PINI)
    if TINI == 0.0:
            Temps     = TP[0]
            Ti        = list(map(float, Temps.tolist()))
            Temps     = Ti 
            ResTime0  = TP[2]
            ResTimei  = [list(map(float, ResTime0))[-1]]
            resTime = ResTimei
    elif TINI >= 1.0:
            Temps    = [ReactorsDict["T_INITIAL"]]; 
            Temps    = list(map(float, Temps)); print(Temps)
            ResTime0 = TP[0]
            ResTimei = list(map(float, ResTime0))
            resTime  = ResTimei; print(ResTimei)
#   except:
#       print("No valid Initial Temperature in Input File...")
#       sys.exit(1)

    if PINI == 0.0:
            press     = TP[1]
            Po        = press
            Pi        = [list(map(float, Po))[0]]
            pressures = Pi
    elif PINI >= 1.0:
            press     = [ReactorsDict["P_INITIAL"]]
            pressures = list(map(float, press)); print(pressures)

#   except:
#       print("No valid Initial Pressure in Input File...")
#       sys.exit(1)
    print("     ","**"*43)
    print("\t> Running case:\n")
    print("\t- List of temperatures:")
    print("\t ",*Temps," / K")
    print("\t- Initial pressure:")
    print("\t ",round(pressures[0],2)," / atm")
    print("\t- Residence time:")
    print("\t ", round(ResTimei[0],2)," / s")
    print("     ","**"*43)


    for m in range(len(ResTimei)):
        if ResTimei[m] == 0:
            ResTimei[m] = 0.00001
        else:
            pass
    Species = list(DataEspeciesJSR)#.split('\\t')
    SP = len(Species); #print(SP)
    #print(Species)
    Prods = [];
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


    #cwd = os.getcwd()
    filenames = Mech


    phi=[ReactorsDict["PHI"]]
    conditions={ReactorsDict["REACTANTS"]}

    speciesToPlot= Prods
    observables=Prods
    volume=1e-6
    sens=sa  #Edit to 1 to run sensitivities
    ###########################################################################
    ###########################################################################
    
    resultsJSR = JSR(filenames,Temps,pressures,resTime,conditions,sens,volume,observables)
    #(resJSR) = Solutions(resultsJSR,
    #                   speciesToPlot,
    #                   conditions,
    #                   pressures,
    #                   phi,
    #                   FinalPath)
    #if sens==1:
    #print("FinalPath: ", FinalPath)
    #print("Data: ")
    #print(observables)
    #print(TP[3:-1])
    DFn       = pd.DataFrame()
    ExpTempos = []
    for u in range(len(observables)):
        dfx = pd.DataFrame({observables[u]:TP[u+3]})
        DFn  = pd.concat([DFn,dfx], axis=1)        
    ExpTempos.append(TP[0])        
    #print(ExpTempos)
    #print(ExpTempos[0])
    a, b, c = plotSpecies(resultsJSR,speciesToPlot,conditions,pressures,phi,FinalPath,NamesDict["InpFileName"],NamesDict["AUTHOR"],DFn,ExpTempos,filenames)
    #else:
    #    pass
    #if sens==1:
    #plotsens(a,observables,conditions,filenames,Temps,4,FinalPath)
    #
    SimTime = time.time() - start
    #resultsJSR.
    #print(b[observables[0]])
    # PLOTTING:
    #for w in range(len(observables)):
    #    plt.rcParams["font.weight"] = "bold"
    #    plt.tick_params(axis= 'both', direction='out', length=4, width=2, labelsize=10)
    #    plt.rcParams["axes.labelweight"] = "bold"
    #    plt.rcParams["axes.linewidth"] = "3"
    #    #plt.scatter(ExpTempos[0],DFn[observables[w]], marker="s", label=str(observables[w].translate(subscript)), lw=4,color='black')
    #    plt.plot(ExpTempos[0],b[observables[w]], label=str(observables[w].translate(subscript)))
    #    #plt.savefig(c + "\\" + NamesDict["InpFileName"] + "_" + observables[w] +'.png',dpi=800, bbox_inches="tight" )
    #    Y  = pd.to_numeric(DFn[observables[w]]) 
    #    YY = (Y).sort_values(ascending=False)
    #    plt.scatter(ExpTempos[0],YY, marker="s", label=str(observables[w].translate(subscript)), lw=4,color='black')
    #    plt.legend()
    #    plt.show()
    #    plt.close()
    #print(b)
    print(">It took {0:0.1f} seconds".format(SimTime))
    return (resultsJSR)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
       (NamesDict, ReactorsDict, UnitsDict,DataEspeciesJSR,dataJSR,FinalPath) = InputFileReaderJSR(Input, "C:\\DKM\\OUTPUTDATA\\");
       JSRCanteraSimulator(NamesDict, ReactorsDict, UnitsDict, DataEspeciesJSR, dataJSR, FinalPath, m, 1)
#============================================================================================
# Eppur si muove! - Galileo Galilei (1564-1642)
#============================================================================================
