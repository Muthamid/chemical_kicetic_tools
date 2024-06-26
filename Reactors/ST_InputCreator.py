# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:29:46 2018
@author: Sergio Martinez
This code snippet extracts experimental data from Excel files and
creates a Shick Tube (ST) input file for running simulations
using the Cantera toolkit.
"""
#===========================================
import os
import glob
import pandas as pd
#============================================

path = os.getcwd(); 

#############################################
# EXTRACTION DATA FROM RCM Excel sheets
#############################################

ExcelName = glob.glob(path+"\\*.xlsx")
if len(ExcelName) == 0:
	print("\t> No excel file found!\n\t> Impossible to create STinputfile\n\t> *.ST_IDT")
else:
	for r in ExcelName:
		ExcelPhi = pd.read_excel(r,sheet_name="Mixture")
	
		AuthorKeyword = "Experimentalist(s):"
		ExcelAuthor = pd.read_excel(r,sheet_name="Read_Me")    
		AUTHOR_INDEX = list(ExcelAuthor).index("Edition:")
		AUTH_EXP     = (ExcelAuthor.iloc[3,1])
	
		ExcelFile = pd.read_excel(r,sheet_name="Exp_Data")    
	
		COL_INDEX      = list(ExcelFile).index("Status")
		COMMENTS_INDEX = list(ExcelFile).index("Comments")
		IDT1st_INDEX   = list(ExcelFile).index("1st Stage_IDT (us)")
		IDT2nd_INDEX   = list(ExcelFile).index("2nd Stage_IDT (us)")
		IDTCRI_INDEX   = list(ExcelFile).index("IDT Criteria")
		TI_INDEX       = list(ExcelFile).index("T5 (K)")
		PI_INDEX       = list(ExcelFile).index("P5 (bar)")
	
		IDT2 = [];IDT1 = []; TII = []; PII=[]; 
		TC = []; PC =[]; CRITERIA = [];
		for u in range(len(ExcelFile)):
			IDT1_EXP     = (ExcelFile.iloc[u,IDT1st_INDEX])
			IDT2_EXP     = (ExcelFile.iloc[u,IDT2nd_INDEX])
			IDT_CRITERIA = (ExcelFile.iloc[u,IDTCRI_INDEX])
			TI_EXP       = (ExcelFile.iloc[u,TI_INDEX])
			PI_EXP       = (ExcelFile.iloc[u,PI_INDEX])
	
			if ExcelFile.iloc[u,COL_INDEX] == 'Accepted':
				IDT1.append(float(IDT1_EXP));
				IDT2.append(float(IDT2_EXP));
				if   "Kistler" in IDT_CRITERIA:
					CRITERIA.append("pressure")
				elif "PMT-CH*" in IDT_CRITERIA:
					CRITERIA.append("ch*")
				elif "PDA-CH*" in IDT_CRITERIA:
					CRITERIA.append("ch*")
				elif "PMT-OH*" in IDT_CRITERIA:
					CRITERIA.append("oh*")
				elif "PDA-OH*" in IDT_CRITERIA:
					CRITERIA.append("oh*")
				else:
					CRITERIA.append(str(IDT_CRITERIA))
				TII.append(float(TI_EXP)); 
				PII.append(float(PI_EXP))

		#----------------------------------
		# BUILDIING THE OUTPUTFILE FORMAT
		#----------------------------------
		specie = []; concen = [];  together = [];
		ExcelMix  = pd.read_excel(r,sheet_name="Mixture")
		Phi_INDEX = list(ExcelMix).index("Phi")
		PHI       = ExcelMix.iloc[0,Phi_INDEX]; 
		P_INDEX   = list(ExcelMix).index("P / bar")
		P         = ExcelMix.iloc[0,P_INDEX]; 
		NameIndex = list(ExcelMix).index("Name")
		ConcIndex = list(ExcelMix).index("[C]")
	
		for y in range(len(ExcelMix)):
			IdtExp     = (ExcelMix.iloc[y,NameIndex])
			ValExp     = (ExcelMix.iloc[y,ConcIndex])
			if "nan" in str(IdtExp) and "nan" in str(ValExp):
				pass
			else:
				specie.append(IdtExp)
				concen.append(ValExp)
				together.append(str(IdtExp)+' '+str(ValExp))
			
		if round(concen[0]) == 0:
			f = open(path+'\\'+str(P)+"_BAR"+'_'+str(round(concen[0],3))+'_'+str(specie[0])+'_PHI_'+str(PHI)+"_"+ ".ST_IDT", "w")
		else:
			f = open(path+'\\'+str(P)+"_BAR"+'_'+str(round(concen[0]))+'_'+str(specie[0])+'_PHI_'+str(PHI)+"_"+ ".ST_IDT", "w")
		f.write("!\AUTHOR: "+str(AUTH_EXP)+"\n")
		f.write("!\JOURNAL: UNPUBLISHED\n")
		f.write("!\PDF: NON\n")
		f.write("!\REACTOR: ST\n")
		f.write("!\DATA: IDT\n")
		f.write("!\PRESSURE: BAR\n")
		f.write("!\TEMPERATURE: K\n")
		f.write("!\TIME: MUSEC\n")
		for m in range(len(specie)):
			f.write('!\\REAC:\t'+str(specie[m])+'\t'+str(concen[m])+'\n')
	
		f.write("!\PHI: "+str(PHI)+"\n")
		f.write("!\IDT_INDICATOR: "+str(CRITERIA[0]).upper()+"_MRC\n")
		f.write("!\IDT_CORRELATION: "+str(specie[0])+"\n")
		if (IDT1[0]) == 0:
			Props = ['p/bar', 'T/K', 'IDT']
			f.write('\t'.join(Props)+'\n')
			for n in range(len(PII)):
				f.write(str(PII[n])+'\t'+str(TII[n])+'\t'+str(IDT2[n])+'\n')
		else:
			Props = ['p/bar','T/K','IDT1','IDT2']
			f.write('\t'.join(Props)+'\n')
			for n in range(len(PII)):
				f.write(str(PII[n])+'\t'+str(TII[n])+'\t'+str(IDT1[n])+'\t'+str(IDT2[n])+'\n')
		f.close()
		print("\t> File name:")
		print("\t. "+str('F'+str(PHI)+"_"+str((P))+"BAR.ST_IDT"))
#==============
# END OF SCRIPT
#==============
