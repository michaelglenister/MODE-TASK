#Calculates and Returns Diagonals of Correlated Matrix for a given set of modes
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import numpy as np
import os, sys

#Lets say that the user has performed NMA on two coarse grained models of the same protein and now wants to compare to see if the additional coarse graining decreased the accuracy. If we obtain the same mean square flucations for each residue then in each model then we can say that the results are comparable regardless of the coarse graining level. But obviously must compare only the residues that are commmon in each model. hence we specify commonResidues (There is a script that works these out as well so we must just link the output of that) here i have specified by chain 

CommonResidues = {'A': [3, 53, 59, 76, 111, 122, 142, 240, 257], 'C': [2, 7, 32, 71, 146], 'B': [52, 56, 70, 142, 161], 'D': [26]} #residues common to multiple pdb files 
interfaceIndex = []
f = open('3VBSProtomer.pdb','r')# The pdb file on which NMA analysis was performed, this script handles on model at a time but we can change this 
NMA = f.readlines()
f.close()
count =0
for line in NMA:
	if line.startswith("ATOM"):
		info = line.split()
		atype = info[2].strip()
		resType = info[3].strip()
		chain = info[4].strip()
		res = int(info[5].strip())
		if (atype == "CB" or (atype=="CA" and resType =="GLY")):
			if res in CommonResidues[chain]: 
				interfaceIndex.append(count) #gets the index of the common atoms, as they would appear in the output W, U and VT matrix. As these matrices contain info on all atoms
			count+=1		
				

Virus = "3VBS_Protomer"
#Specify modes
TotalModes = 2526 #number of residues in protein *3
FirstMode = 2519  #input user
LastMode = 2519 #input user

#Specify Residue Indexes
FirstRes = 1
LastRes = 842
ResRange = range(FirstRes-1,LastRes)

ModeRange = range(FirstMode,LastMode+1)
print ModeRange

fw = open('3VBSProtomerW.txt','r') #W matrix input file that was output from C++ Scripts
EigenValues = fw.readlines()
fw.close()

#Create A Full W Inverse Matrix (This is if we want correlation averaged over all modes)
WInv = np.zeros((TotalModes,TotalModes))
for i in range(TotalModes):
	if i<TotalModes-6:
                WInv[i,i] = 1/(float(EigenValues[i].split()[1].strip())) 
print "W in "


#Create Filtered W Inverse Matrix (This is if we want correlation for a specific mode)
WF = np.zeros((TotalModes,TotalModes))
for i in ModeRange:
	WF[i,i] = 1/(float(EigenValues[i].split()[1].strip()))
print "WF in "

#Read In U and VT full Matrix as U is the transpose of VT I only read in VT and create U from the VT matrix info. So we can exlcude U output from C++ script for faster analysis 
fvt = open('3VBSProtomerVT.txt','r')
EigenVectors = fvt.readlines()
fvt.close()
print "U and VT file in "

VT = np.zeros((TotalModes,TotalModes))
U = np.zeros((TotalModes,TotalModes))

for i in range(TotalModes):
	vectors = EigenVectors[i].split()
	for j in range(TotalModes):
		vector = float(vectors[j].strip())
		VT[i,j] = vector
		U[j,i] = vector

print "U and VT read"
#Calculate Correlation Matrices
#Full C matrix
WVT = np.dot(WInv,VT)
print "Correlations Calculated"
C = np.dot(U,WVT)
print "Correlations Calculated"

#Mode Specific C Matrix
WVTm = np.dot(WF,VT)
print "Correlations Calculated"
CM = np.dot(U,WVTm)

print "Correlations Calculated"


#Calculate Trace of the Correlation Matrices
TraceC = np.zeros((TotalModes/3,TotalModes/3))
TraceCM = np.zeros((TotalModes/3,TotalModes/3))

for i in range(0,TotalModes,3):
	for j in range(0,TotalModes,3):
		trace = 0
		traceM = 0
		for k in range(3):
			trace = trace+C[i+k,j+k]
			traceM = traceM+CM[i+k,j+k]
		TraceC[i/3,j/3]=trace
		TraceCM[i/3,j/3]=traceM

#Print the diagonal values per residue


w = open(Virus+str(FirstMode)+"BetaValues.txt",'w')
w.write("Full Correaltion\n")
w.write("Res\tBetaValue\n")
for i in ResRange:
	w.write(str(i+1)+"\t"+str(TraceC[i,i])+"\n")
w.write("Common Residues")
for i in interfaceIndex:
	w.write(str(i+1)+"\t"+str(TraceC[i,i])+"\n")


w.write("\nFiltered Correaltion\n")
w.write("Res\tBetaValue\n")
for i in ResRange:
	w.write(str(i+1)+"\t"+str(TraceCM[i,i])+"\n")
w.write("Common Residues")
for i in interfaceIndex:
	w.write(str(i+1)+"\t"+str(TraceCM[i,i])+"\n")
w.close()
















