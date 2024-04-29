# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:29:46 2018
@author: Sergio Martinez
This code snippet extracts the experimental data from Excel files and
creater a Rapid Compression Machine (RCM) input file to run simulations
using the Cantera toolkit.
"""
#===========================================
import os
import sys
import glob
import numpy as np
import pandas as pd
import cantera as ct
#============================================

path = os.getcwd(); 
form = '\\*.VPRO'
#--------------------------------------------
# FOLDERS WITH VPRO FILES
#--------------------------------------------

VFNames = []; Ti = []; Pi = [];
FOLDER = glob.glob(path+form);
if len(FOLDER) == 0:
	New = glob.glob(path+'\\*.inp_r');
	for v in range(len(New)):
		newname = New[v].replace('.inp_r', '.VPRO');
		output = os.rename(New[v], newname)
	vpros = glob.glob(path+'\\*.VPRO')
	if len(vpros) == 0:
		for v in range(len(New)):
			newname = New[v].replace('.INP_r', '.VPRO');
			output = os.rename(New[v], newname)
	vpros = glob.glob(path+'\\*.VPRO')
	if len(vpros) == 0:
		for v in range(len(New)):
			newname = New[v].replace('.INP_R', '.VPRO');
			output = os.rename(New[v], newname)
	vpros = glob.glob(path+'\\*.VPRO')
	if len(vpros) == 0:
		for v in range(len(New)):
			newname = New[v].replace('.inp_R', '.VPRO');
			output = os.rename(New[v], newname)
else:
	pass

form = '\\*.VPRO';
FOLDER = glob.glob(path+form);
FileReac = pd.read_csv(FOLDER[0],names=['col']);
Foldername = (FOLDER[0].split('\\'))[-2]
REACS = ((FileReac[FileReac.col.str.contains("REAC")==True]).col.str.split());

#############################################
# EXTRACTION DATA FROM RCM Excel sheets
#############################################
ExcelName = glob.glob(path+"\\*.xlsx")
ExcelPhi = pd.read_excel(ExcelName[0],sheet_name="Mixture")

AuthorKeyword = "Experimentalist(s):"
ExcelAuthor = pd.read_excel(ExcelName[0],sheet_name="Read_Me")    
AUTHOR_INDEX = list(ExcelAuthor).index("Edition:")
AUTH_EXP     = (ExcelAuthor.iloc[3,1])

ExcelFile = pd.read_excel(ExcelName[0],sheet_name="Exp_Data")    

COL_INDEX      = list(ExcelFile).index("Status")
COMMENTS_INDEX = list(ExcelFile).index("Comments")
IDT2nd_INDEX   = list(ExcelFile).index("2nd Stage_IDT (ms)")
IDT1st_INDEX   = list(ExcelFile).index("1st Stage_IDT (ms)")
IDTCRI_INDEX   = list(ExcelFile).index("IDT Criteria")
TEMP_INDEX     = list(ExcelFile).index("Tc (K)")
PRESS_INDEX    = list(ExcelFile).index("Pc (bar) Kistler + Pi")
TI_INDEX       = list(ExcelFile).index("Ti (K)")
PI_INDEX       = list(ExcelFile).index("Pi (bar)")

IDT2 = [];IDT1 = []; TII = []; PII=[]; 
TC = []; PC =[]; CRITERIA = [];
for u in range(len(ExcelFile)):
	IDT2_EXP     = (ExcelFile.iloc[u,IDT2nd_INDEX])
	IDT1_EXP     = (ExcelFile.iloc[u,IDT1st_INDEX])
	IDT_CRITERIA = (ExcelFile.iloc[u,IDTCRI_INDEX])
	TEMP_EXP     = (ExcelFile.iloc[u,TEMP_INDEX])
	PRESS_EXP    = (ExcelFile.iloc[u,PRESS_INDEX])	
	TI_EXP       = (ExcelFile.iloc[u,TI_INDEX])
	PI_EXP       = (ExcelFile.iloc[u,PI_INDEX])
    
	if ExcelFile.iloc[u,COL_INDEX] == 'Accepted':
		TC.append(float(TEMP_EXP)); 
		PC.append(float(PRESS_EXP))
		IDT2.append(float(IDT2_EXP));
		IDT1.append(float(IDT1_EXP));
		CRITERIA.append(str(IDT_CRITERIA))
	#if ExcelFile.iloc[u,COMMENTS_INDEX] == 'NR':
		TII.append(float(TI_EXP)); 
		PII.append(float(PI_EXP))



#============================================
# EXTRACTION DATA FROM RCM VPRO FILES
#============================================
print("\t> Using the next volume profile files (VPRO): \n\t")
for i in range(len(FOLDER)):
	VproFileNames = (FOLDER[i].split('\\'))[-1]; print("\t. "+VproFileNames)
	VFNames.append(VproFileNames);
	FileRCM = pd.read_csv(FOLDER[i],names=['col']);
	TI = ((FileRCM[FileRCM.col.str.contains("TEMP")==True]).col.str.split());
	TI = float(TI.iloc[0][1])
	Ti.append(TI);# print(TI)
	PI = ((FileRCM[FileRCM.col.str.contains("PRES")==True]).col.str.split());
	PI = (float(PI.iloc[0][1]))
	Pi.append(PI);# print(PI)
print("\n")

#----------------------------------
# BUILDIING THE OUTPUTFILE FORMAT
#----------------------------------
#print("\t> Using the next volume profile files (VPRO): \n\t"+str(len(VFNames)))
#ry:
#	FoldNamecrack = (Foldername.split('_'))[0]
#	PHI = (FoldNamecrack.split('F'))[-1] ; Phis = float(PHI)/100
#except:
PHI = ExcelPhi.iloc[4,3]; print(PHI)
Phis = str(PHI).split(" ")[2]
#print(Phis)
f = open(path+'\\F'+str(Phis)+"_"+str(int(PC[0]))+"BAR.RCM_IDT", "w")
f.write("!\AUTHOR: "+str(AUTH_EXP)+"\n")
f.write("!\JOURNAL: N/A\n")
f.write("!\PDF: N/A\n")
f.write("!\REACTOR: RCM\n")
f.write("!\DATA: IDT\n")
f.write("!\PRESSURE: BAR\n")
f.write("!\TEMPERATURE: K\n")
f.write("!\TIME: MSEC\n")
f.write("!\PHI: "+str(Phis)+"\n")
for m in range(len(REACS)):
    f.write('!\\REAC:\t'+str((REACS.iloc[m][1]+'\t'+REACS.iloc[m][2]))+'\n')
f.write("!\IDT_INDICATOR: "+str(CRITERIA[0]).upper()+"_MRC\n")
f.write("!\IDT_CORRELATION: "+str(REACS.iloc[0][1])+"\n")
f.write("!\COMMENT: ONLY HAVE REPRESENTATIVE NON-REACTIVE, 0-VALUE IDT CONTIDIONS WILL NOT BE PLOTTED\n")
Props = ['TI','PI','TC','PC','IDT_STG1','IDT_STG2','FACILITY_EFFECT']
f.write('\t'.join(Props)+'\n')
Tc = []; Pc =[]; STG1 = []; STG2 = [];
for k in range(len(TC)):
    Tc.append(round(TC[k],2))
    Pc.append(round(PC[k],2))
    STG1.append(IDT1[k])
    STG2.append(IDT2[k])
#print(len(TC))
for p in range(len(TC)):
	cwd = os.getcwd()
	vpro = glob.glob(cwd+"\\*"+str(int((TII[p])))+"*.VPRO"); print(PII[p])
	if len(vpro) > 0:
		vpro = str(vpro).split("\\")[-1]; print(vpro)
		vpro = vpro.split("'")[0]; print(vpro)
		f.write((str(round(TII[p],3)))+'\t'+str(round(PII[p],3))+'\t'+str(Tc[p])+'\t'+str(Pc[p])+'\t'+str(STG1[p])+'\t'+str(STG2[p])+"\t"+str(vpro)+"\n")
	else:
		vpro = glob.glob(cwd+"\\*"+str(int((PII[p])))+"*.VPRO"); print(PII[p])
		vpro = str(vpro).split("\\")[-1]; print(vpro)
		vpro = vpro.split("'")[0]; print(vpro)
		f.write((str(round(TII[p],3)))+'\t'+str(round(PII[p],3))+'\t'+str(Tc[p])+'\t'+str(Pc[p])+'\t'+str(STG1[p])+'\t'+str(STG2[p])+"\t"+str(vpro)+"\n")
		
f.close()
print("\t> File name:")
print("\t. "+str('F'+str(Phis)+"_"+str(int(PC[0]))+"BAR.RCM_IDT"))
#==============
# END OF SCRIPT
#==============
