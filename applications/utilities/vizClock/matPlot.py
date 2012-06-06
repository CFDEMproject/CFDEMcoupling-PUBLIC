#!/usr/bin/env python
import csv, sys
import numpy as np
import matplotlib.pyplot as plt

# Open the data
datafile = "timeEvalFull.txt"
f = open(datafile, 'r')
reader = csv.reader(f, dialect='excel-tab')
reader.next()

header = []
identifier = []
deltaT = []
maxdeltaT = []
nOfRuns = []
level = []
parentNr = []
parentName = []


i = 0
for row in reader:
	if i == 0:
		for column in row:
			header.append(column)
		print header
	else:
		identifier.append(row[0])
		deltaT.append(float(row[1]))
		maxdeltaT.append(float(row[2]))
		nOfRuns.append(int(row[3]))
		level.append(int(row[4]))
		parentNr.append(int(row[5]))
		parentName.append(row[6])
	i+=1

print identifier
print deltaT
print maxdeltaT
print nOfRuns
print level
print parentNr
print parentName


bottom = []
brotherheight = []
for i in range(len(identifier)):
	bottom.append(0)
	brotherheight.append(0)

for i in range(len(identifier)):
	if i != 0:
		if level[i]<level[i-1]:
			brotherheight[level[i-1]]=0
	if parentNr[i] != -1:
		bottom[i] = bottom[parentNr[i]]
	bottom[i] += brotherheight[level[i]]
	brotherheight[level[i]] += deltaT[i]

for i in range(len(identifier)):
	plt.bar(level[i],deltaT[i],width = 0.2, bottom=bottom[i])
	plt.text(level[i]+0.22,bottom[i]+deltaT[i]/2,identifier[i]+" "+str(nOfRuns[i])+"x")
plt.xlabel('run level')
plt.ylabel('CPU time in s')
plt.title('time measurement')
plt.show()


