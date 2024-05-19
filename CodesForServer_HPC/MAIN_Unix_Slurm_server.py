#########################################################################
# " Eppur si muove, Sinead na Paor - S.M. (2020-3986) "   #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Import Modules
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import sys
import os
import stat
import subprocess
from subprocess import call
#from itertools import izip_longest
#from DosTitles.DosTitles import PrintDosHeader
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#print"Edited paths and ICHEC project by VPatel on 01/08/2021 "
#print"\n"
#print"/|\\"*25
#print" All cleaned!                                   "
#print" Thanks for use this Code                       "
#print" HOPE YOU ENJOY IT!                             "
#print"                                                "
#print" Guv + A*(guv) = [8*(pi)*G]/(Tuv)               "
#print" Hope all your problems stay beyond             "
#print" the boundary of the event horizon,             "
#print" so they won't affect you any more...           "
#print" Gral. Rel. Eq. Albert Einstein                 "
#print"                                                "
#print" @JMSPSoft@ ^ >.< ^ @SineadnaPoir Subroutine@   "
#print"               U                                "
#print"             |___|                              "
#print"              \_/                               "
#print" S.M. 16/01/2020                                "
#print"/|\\"*25
################################
## 15/01/2020  Sergio Martinez##
################################"""
"""
example(s):

	f = open("jobST.sh","w+")
	f.write("#!/bin/sh\n")
	f.write("#SBATCH -t 01:00:00\n")
	f.write("#SBATCH -p DevQ\n")
	f.write("#SBATCH -N 1\n")
	f.write("#SBATCH -A ngche114c\n")
	f.write("#SBATCH --job-name=ST\n")
	f.write("\n")
	f.write("cd $SLURM_SUBMIT_DIR\n")
	f.write("\n")
	f.write("module load taskfarm\n")
	f.write("module load conda/2\n")
	f.write("source activate val\n")
	f.write("\n")
	f.write("time taskfarm tasksST\n")
	f.close()

    f = open("jobRCM.sh","w+")
    f.write("#!/bin/sh\n")
    f.write("#SBATCH -t 05:00:00\n")
    f.write("#SBATCH -p ProdQ\n")
    f.write("#SBATCH -N 1\n")
    f.write("#SBATCH -A ngche114c\n")
    f.write("#SBATCH --job-name=RCM\n")
    f.write("\n")
    f.write("cd $SLURM_SUBMIT_DIR\n")
    f.write("\n")
    f.write("module load taskfarm\n")
    f.write("module load conda/2\n")
    f.write("source activate val\n")
    f.write("\n")
    f.write("time taskfarm tasksRCM\n")
    f.close()

"""
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Reads in the bin control file that contains paths to databases/folders that
# are regularly used -- takes as input the name of the BinCon file
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def ReadBinControlFile(file_path,inp_file_bin):
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Read in the BinCon file
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    PrintDosHeader("", "=", 90)
    print(">Reading BinCon File:", file_path)
    PrintDosHeader("", "=", 90)
    dict = {
        "ExpsPath": "/ichec/work/ngche114c/vpatel/DKM/DATA/",
        "ModelsPath": "/ichec/work/ngche114c/vpatel/DKM/MODELS/",
        "OutputPath": "/ichec/work/ngche114c/vpatel/DKM/OUTPUTDATA/",
        "ExpDB": inp_file_bin,
        "ModelDB": inp_file_bin
        }
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Open the file, read it in, and provide
    # a dictionary containing keys that provide
    # paths to certain bins
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    return(dict)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Gets a list of subdirectories and files within the path provided
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def GetDirFileLists(base_path):
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Walk the dir and return lists of files, dirs etc.
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    returnval = "Yes"
    DirList = []
    FileList = []
    for root, dirs, files in os.walk(base_path, topdown=True):
        DirList.append(root)
        for f in files:
            FileList.append(os.path.join(root, f))
    return(DirList, FileList)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CWD           = "/ichec/work/ngche114c/vpatel/DKM/DATA/"  
if (len(sys.argv)) == 5:
	ListOfSpecies = str(sys.argv[2]).rsplit(',')  
	ST   = []
	RCM  = []
	LBV  = []
	JSR  = []
	FR   = []
	BIN  = []
	DIRS = []; FILES = []; WPATH = []
	for k in ListOfSpecies:
		WPath      = CWD+k; WPATH.append(WPath)
		Dirs,Files = GetDirFileLists(str(WPath))
		for j in range(len(Dirs)):
			DIRS.append(Dirs[j]);
		for j in range(len(Files)):
			FILES.append(Files[j])
	# Reading all the data in directories/subdirs inside 
	counter = 1
	for i in FILES:
        	if  i.endswith('.ST_IDT'):
                	ST.append(i)
        	elif i.endswith('.RCM_IDT'):
                	RCM.append(i)
        	elif i.endswith('.FPF_EXP'):
                	LBV.append(i)
        	elif i.endswith('.JSR_EXP'):
                	JSR.append(i)
        	elif i.endswith('.FR_EXP'):
                	FR.append(i)
        	else:
                	BIN.append(i)
	Dirs = DIRS; Files = FILES; WPath = WPATH
else:
	Dirs,Files = GetDirFileLists(CWD)
	ST  = []
	RCM = []
	LBV = []
	JSR = []
	FR  = []
	BIN = []
	counter = 1
	for i in Files:
		if  i.endswith('.ST_IDT'):
			ST.append(i)
		elif i.endswith('.RCM_IDT'):
			RCM.append(i)
		elif i.endswith('.FPF_EXP'):
			LBV.append(i)
		elif i.endswith('.JSR_EXP'):
			JSR.append(i)
		elif i.endswith('.FR_EXP'):
			FR.append(i)
		else:
			BIN.append(i)
	WPath = [CWD]
print"-.-"*20
print"\tRunning SdP... "
print"\tICHEC version1.0"
print"\tS.M.-C3-NUIGalway, IE"
print"\t0b11111100100"
print"\t01001010,01001101,01010011,01010000"
print"\t0b11111010101"
print"-.-"*20
print"\n"
print"\t> Working path given:"
for b in WPath:
	print"\t"+str(b)
print"+"*60
print"\t> Directories: ",len(Dirs)
print"\t> Files: ",len(Files)
print"+"*60
print"+"*60
print"\t> ST  files: ",len(ST)
print"\t> RCM files: ",len(RCM)
print"\t> LBV files: ",len(LBV)
print"\t> JSR files: ",len(JSR)
print"\t> FR  files: ",len(FR)
print"\t> Total Reac files: ",len(ST)+len(RCM)+len(LBV)+len(JSR)+len(FR)
print"\t  Others formats : ",len(BIN)
print"\t  Total files scanned: ",len(ST)+len(RCM)+len(LBV)+len(JSR)+len(FR)+len(BIN)
print"+"*60
if (len(sys.argv)) > 1:
 if   str(sys.argv[1]) == "ST":
      # Generating all slurm files to run individual simulations
      # Making a contol bash script file to run all the simulations in one go...
      f = open("jobST_"+str(sys.argv[3])+".sh","w+")
      #f = open("jobST.sh","w+")
      f.write("#!/bin/sh\n")
      f.write("#SBATCH -t 1:00:00\n")
      f.write("#SBATCH -p DevQ\n")
      f.write("#SBATCH -N 1\n")
      f.write("#SBATCH -A ngche114c\n")
      f.write("#SBATCH --job-name=ST\n")
      f.write("\n")
      f.write("cd $SLURM_SUBMIT_DIR\n")
      f.write("\n")
      f.write("module load taskfarm\n")
      f.write("module load conda/2\n")
      f.write("source activate val\n")
      f.write("\n")
      f.write("time taskfarm tasksST_"+str(sys.argv[3]))
      #f.write("time taskfarm tasksST\n")
      f.close()
      # Giving privilegies to st.sh file, to execute from python
      # Making a contol bash script file to run all the simulations in one go...
      f = open("tasksST_"+str(sys.argv[3]),"w+")
      #f = open("tasksST","w+")
      f.close()
      for i in ST:
          subprocess.call("./probST_"+str(sys.argv[3])+".sh "+i+" "+str(sys.argv[4]),shell=True)
          #subprocess.call("./probST.sh "+i+" "+str(sys.argv[4]),shell=True)
      ### Here add the calling to BASH to run all simulations prepared before
      f = open("callST.sh","w+")
      f.write("#!/bin/sh\n")
      f.write("sbatch jobST_"+str(sys.argv[3])+".sh")
      #f.write("sbatch jobST.sh")
      f.close()
      other = os.stat("callST.sh")
      os.chmod("callST.sh",other.st_mode | stat.S_IEXEC)
      subprocess.call("./callST.sh")
      os.remove("callST.sh")
 #####
 # RCM
 #####
 elif str(sys.argv[1]) == "RCM":
      # Generating all slurm files to run individual simulations
      # Making a contol bash script file to run all the simulations in one go... ProdQ/DevQ
      f = open("jobRCM_"+str(sys.argv[3])+".sh","w+")
      f.write("#!/bin/sh\n")
      f.write("#SBATCH -t 8:00:00\n")
      f.write("#SBATCH -p ProdQ\n")
      f.write("#SBATCH -N 1\n")
      f.write("#SBATCH -A ngche114c\n")
      f.write("#SBATCH --job-name=RCM\n")
      f.write("\n")
      f.write("cd $SLURM_SUBMIT_DIR\n")
      f.write("\n")
      f.write("module load taskfarm\n")
      f.write("module load conda/2\n")
      f.write("source activate val\n")
      f.write("\n")
      f.write("time taskfarm tasksRCM_"+str(sys.argv[3]))
      f.close()
      #print(str(sys.argv[3]))
      # Giving privilegies to st.sh file, to execute from python
      # Making a contol bash script file to run all the simulations in one go...
      f = open("tasksRCM_"+str(sys.argv[3]),"w+")
      f.close()
      ### Reading all the data in directories/subdirs inside "/ichec/work/ngchec079c/sergito/DKM/DATA"
      for i in RCM:
          subprocess.call("./probRCM_"+str(sys.argv[3])+".sh "+i+" "+str(sys.argv[4]),shell=True)
      ### Here add the calling to BASH to run all simulations prepared before
      f = open("callRCM.sh","w+")
      f.write("#!/bin/sh\n")
      f.write("sbatch jobRCM_"+str(sys.argv[3])+".sh")
      f.close()
      other = os.stat("callRCM.sh")
      os.chmod("callRCM.sh",other.st_mode | stat.S_IEXEC)
      subprocess.call("./callRCM.sh")
      os.remove("callRCM.sh")
 # ###
 # LBV
 # ###
 elif str(sys.argv[1]) == "LBV":
      # Generating all slurm files to run individual simulations
      # Making a contol bash script file to run all the simulations in one go... ProdQ/DevQ
      f = open("jobLBV_"+str(sys.argv[3])+".sh","w+")
      #f = open("jobLBV.sh","w+")
      f.write("#!/bin/sh\n")
      f.write("#SBATCH -t 72:00:00\n")
      f.write("#SBATCH -p ProdQ\n")
      f.write("#SBATCH -N 1\n")
      f.write("#SBATCH -A ngche114c\n")
      f.write("#SBATCH --job-name=LBV\n")
      f.write("\n")
      f.write("cd $SLURM_SUBMIT_DIR\n")
      f.write("\n")
      f.write("module load taskfarm\n")
      f.write("module load conda/2\n")
      f.write("source activate val\n")
      f.write("\n")
      f.write("time taskfarm tasksLBV_"+str(sys.argv[3]))
      #f.write("time taskfarm tasksLBV\n")
      f.close()
      # Giving privilegies to st.sh file, to execute from python
      # Making a contol bash script file to run all the simulations in one go...
      f = open("tasksLBV_"+str(sys.argv[3]),"w+")
      #f = open("tasksLBV","w+")
      f.close()
      ### Reading all the data in directories/subdirs inside "/ichec/work/ngchec101c/vpatel/DKM/DATA"
      for i in LBV:
          subprocess.call("./probLBV_"+str(sys.argv[3])+".sh "+i+" "+str(sys.argv[4]), shell=True)
          #subprocess.call("./probLBV.sh "+i,shell=True)
      ### Here add the calling to BASH to run all simulations prepared before
      f = open("callLBV.sh","w+")
      f.write("#!/bin/sh\n")
      f.write("sbatch jobLBV_"+str(sys.argv[3])+".sh")
      #f.write("sbatch jobLBV.sh")
      f.close()
      other = os.stat("callLBV.sh")
      os.chmod("callLBV.sh",other.st_mode | stat.S_IEXEC)
      subprocess.call("./callLBV.sh")
      os.remove("callLBV.sh")
 #
 #
 #
 else:
       print("#"*75)
       print("\t SERGITO WARNING!!! PAY ATENTION!!!")
       print("#"*75)
       print("\t> Nothing to do here...")
       print("\t> You need to specify the reactor")
       print("\t> Ex. > python DirTools.py ST")
       #print(
       print("\t> You are welcome, now run it again :)")
       print("-"*75)
else:
	print("#"*75)
	print("\t SERGITO WARNING!!! PAY ATENTION!!!")
	print("#"*75)
	print("\t> Nothing to do here...")
	print("\t> You need to specify the reactor")
	print("\t> Ex. > python DirTools.py ST")
#############################################
## ALL YOU NEED IS LOVE! & Science :)  S.M. #
#############################################
