#!/bin


import numpy as np
from scipy.io import savemat
import matplotlib.pyplot as plt
import shutil
import os
import csv
import sys
import re

import argparse
parser = argparse.ArgumentParser(prog="postProcess.py",
								description="Tornado Simulator OpenFOAM data pot-processing")
parser.add_argument('-r','--resid',action='store_true',
					help='Flag to plot solution residuals')
parser.add_argument('-i','--inlet',action='store_true',
					help='Flag to plot mean velocity history over inlet area')
parser.add_argument('-o','--outflow',action='store_true',
					help='Flag to plot mean velocity history in outflow duct')
parser.add_argument('-p','--profile',action='store_true',
					help='Flag to plot mean vortex velocity profile')
parser.add_argument('-c','--calcSwirlRatio',action='store_true',
					help='Flag to calculate mean swirl ratio')
parser.add_argument('-s',metavar="startTime",type=float,default=0,
					help='Starting Iteration for swirl ratio averaging window')
parser.add_argument('-e',metavar="endTime",type=float,default=-1,
					help='Ending Iteration for swirl ratio averaging window')
parser.add_argument('-nPts',metavar='count',type=int,default=100,
					help='Quantity of points in window to average swirl ratio over (default = 100)')
parser.add_argument('-np',action='store_true',
					help='Flag to disable all plots while debugging')
parser.add_argument('-d',metavar='directory',type=str,default="postProcessing",
					help='postProcessing directory to read from')
parser.add_argument('-t',metavar='time',type=str,default='0',
					help='Start time of OpenFOAM run')
parser.add_argument('--mat',action='store_true',
					help='Flag to save plot data as a .mat file for plotting in MATLAB')
args = parser.parse_args()

def checkArgs():

	if not os.path.isdir(args.d):
		print("Error: Directory does not exist.")
		if args.d == "postProcessing":
			print("    Did you remember to specify a directory?")
		return 0

	return 1




def whichFile(fname):
	# returns the file with the largest file size in case a partial
	#	run was started at a time

	islash = -1
	for i in range(len(fname)):
		j = len(fname)-i-1
		if fname[j] == '/' and islash == -1:
			islash = j
	files = os.listdir(fname[:islash])
	if len(files) == 1:
		return fname
	else:
		fsize = np.zeros(len(files))
		for i in range(len(files)):
			ffile = fname[:islash]+'/'+files[i]
			stat = os.stat(ffile)
			fsize[i] = stat.st_size
		useidx = np.argmax(fsize)
		return fname[:islash]+'/'+files[useidx]


def residualPlot(): # Plots solver U and p residuals
	# fname = "./pp5/solverInfo/0/solverInfoNew.dat"
	fname = "./" +args.d+ "/solverInfo/" + str(args.t) + "/solverInfo.dat"
	fname = whichFile(fname)
	file = open(fname,'r')
	lines = file.readlines()
	iMax = 0
	for line in lines:
		f = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line)
		if len(f) != 0:
			iMax = int(f[0])

	it = np.zeros(iMax)
	Uxr = np.zeros(iMax)
	Uyr = np.zeros(iMax)
	Uzr = np.zeros(iMax)
	pr = np.zeros(iMax)

	i = 0
	for line in lines:
		f = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line)
		if len(f) != 0:
			it[i] = int(f[0])
			Uxr[i] = float(f[1])
			Uyr[i] = float(f[4])
			Uzr[i] = float(f[7])
			pr[i] = float(f[10])
			i += 1

	file.close()

	plt.figure(152)
	plt.yscale('log')
	plt.plot(it,Uxr)
	plt.plot(it,Uyr)
	plt.plot(it,Uzr)
	plt.plot(it,pr)
	plt.grid()
	plt.legend(['U_x','U_y','U_z','p'],loc='upper right')
	plt.xlabel('Iteration')
	plt.ylabel('Initial Residual')


def turn(): # Plots mean velocity in outflow duct turning AD
	fname = "./" +args.d+ "/FOturnMAG/" + str(args.t) + "/volFieldValue.dat"
	fname = whichFile(fname)
	file = open(fname,'r')
	lines = file.readlines()
	for line in lines:
		f = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line)
		if len(f) > 1:
			iMax = int(f[0])

	it = np.zeros(iMax)
	uy = np.zeros(iMax)
	i = 0
	for line in lines:
		f = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line)
		if len(f) > 1:
			it[i] = int(f[0])
			uy[i] = float(f[1])
			i += 1
	print(np.mean(uy[6500:7500]))

	plt.figure(35487)
	plt.clf()
	plt.plot(it,uy)
	plt.grid()
	plt.xlabel("Iteration")
	plt.ylabel("Turning AD Mean Velocity")


def inletGrab(): # Plots mean velocity in main AD inlet
	fname = "./" +args.d+ "/FOinletMAG/" + str(args.t) + "/volFieldValue.dat"
	fname = whichFile(fname)
	file = open(fname,'r')
	lines = file.readlines()
	nlines = 0

	iMin = 0
	for line in lines:
		f = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line)
		if len(f) > 1:
			iMax = int(f[0])
			if f[0] != 0 and iMin == 0:
				iMin = int(f[0])

	it = np.zeros(iMax-iMin+1)
	uy = np.zeros(iMax-iMin+1)
	i = 0
	for line in lines:
		f = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",line)
		if len(f) > 1:
			it[i] = int(f[0])
			uy[i] = float(f[1])
			i += 1

	return [it,uy]

def inletPlot():
	[it,uy] = inletGrab()
	plt.figure(1386)
	plt.clf()
	plt.plot(it,uy)
	plt.grid()
	plt.xlabel("Iteration")
	plt.ylabel("Main AD Mean Velocity")


def getSwirls(fname):
	trt = np.sqrt(2)/2
	data = np.genfromtxt(fname, delimiter=',',skip_header=1)
	r = np.sqrt(np.power(data[:,0],2)+np.power(data[:,2],2))
	r[0:int(len(r)/2)] *= -1
	nx = -data[:,2] / r
	nz = data[:,0] / r
	ur = np.abs(data[:,3]*nx + data[:,5]*nz)

	rcf = 0
	for i in range(len(ur)):
		rcf += ur[i]
	rcf = rcf/len(ur)

	mini = np.argmin(ur)
	rMin = r[mini]
	rMaxL = 0
	rMaxR = 0
	maxL = 0
	maxR = 0

	for i in range(mini):
		if ur[i] > maxL:
			maxL = ur[i]
			rMaxL = np.abs(r[i] - rMin)

	for i in range(mini,len(ur)):
		if ur[i] > maxR:
			maxR = ur[i]
			rMaxR = np.abs(r[i] - rMin)

	maxA = (maxL + maxR) /2
	rMaxA = (rMaxL + rMaxR) / 2

	return [maxA, rMaxA, rcf]

def getVtxFnames():
	# Get vortex files/file names
	files = ["01","02","03","04","11","12","13","14","21","22","23","24","31","32","33","34"]
	nFiles = len(files)
	for i in range(nFiles):
		files[i] = "l" + files[i] + "_U.csv"
	vtxTimes = os.listdir("./" +args.d+ "/FOvortex")
	vtxTime = np.zeros(len(vtxTimes))
	for i in range(len(vtxTimes)):
		vtxTime[i] = float(vtxTimes[i])
	del vtxTimes
	vtxTime = np.sort(vtxTime)
	
	# Shrink window, start/end
	if args.e == -1:
		iEnd = len(vtxTime)
	else:
		foundFlag = 0
		i = len(vtxTime)-1
		while foundFlag == 0:
			if vtxTime[i] > args.e:
				i -= 1
			else:
				iEnd = i + 1
				foundFlag = 1

			if i == 0:
				print("Error: End time not found in time directories")
				break

	if args.s == 0:
		iStart = 0
	else:
		foundFlag = 0
		i = 0
		while foundFlag == 0:
			if vtxTime[i] < args.s:
				i += 1
			else:
				iStart = i
				foundFlag = 1

			if i == len(vtxTime):
				print("Error: End start not found in time directories")
				break

	vtxTime = vtxTime[iStart:iEnd]

	# Shrink window, nPoints
	samples = (len(vtxTime)-1)/args.nPts
	if np.floor(samples) != samples:
		samples = np.floor(samples)
		if samples == 0:
			samples = 1
			
	vtxTime = vtxTime[::int(samples)]

	return [files,vtxTime]


def vp():
	[files,vtxTime] = getVtxFnames()

	profile = np.zeros(250)
	n = 0

	for time in vtxTime:
		if np.floor(time) == time:
			fPath = "./" +args.d+ "/FOvortex/" + str(int(time)) + "/"
		else:
			fPath = "./" +args.d+ "/FOvortex/" + str(time) + "/"

		for i in range(len(files)):
			fname = fPath + files[i]
			trt = np.sqrt(2)/2
			data = np.genfromtxt(fname, delimiter=',',skip_header=1)
			r = np.sqrt(np.power(data[:,0],2)+np.power(data[:,2],2))
			r[0:int(len(r)/2)] *= -1
			nx = -data[:,2] / r
			nz = data[:,0] / r
			ur = np.abs(data[:,3]*nx + data[:,5]*nz)

			if len(ur) == len(profile):
				profile = (profile*n + ur)/(n+1)
				n += 1

	if args.mat:
		if not args.d == "postProcessing":
			rname = 'r_' + args.d
			pname = 'p_' + args.d
			fname = 'vortexProfile_' + args.d + '.mat'
		else:
			rname = 'r'
			pname = 'p'
			fname = 'vortexProfile.mat'
		matdict = {rname: r, pname: profile}
		savemat(fname,matdict)
		print("Saved .mat file with vortex profile.")


	if args.profile:
		plt.figure(6587)
		plt.plot(r[124:],profile[124:])
		plt.xlabel("r (m)")
		plt.ylabel(r'$U_\theta$ (m/s)')
		plt.grid()


def csr():

	[qit,uy] = inletGrab()
	r2in = 0.915**2

	miv = 0 # Mean Inlet Velocity Magnitude
	qiv = 0 # Mean Inlet Mass Flow
	mur = 0 # Mean Utheta Max
	mur_std = 0 # /\ standard deviation
	mrr = 0 # Mean R Utheta Max
	mrr_std = 0 # /\ stdeviation
	msr = 0 # Mean swirl ratio
	msr_std = 0 # /\ stdev

	[files,vtxTime] = getVtxFnames()

	def update(newVal,mean,std,N):
		newMean = ((mean*N) + newVal)/(N+1)
		if N != 0:
			std = np.sqrt(((N-1)*(std*std) + (newVal - newMean)*(newVal - mean))/N)
		return [newMean,std]

	nDirs = len(vtxTime)
	nOverall = 0
	for time in vtxTime:
		nTime = 0
		tur = 0 # Time Utheta Max
		tur_std = 0 # /\ stdeviation
		trr = 0 # Time R Utheta Max
		trr_std = 0 # /\ stdeviation
		tsr = 0 # Mean swirl ratio
		tsr_std = 0 # /\ stdev

		if np.floor(time) == time:
			fPath = "./" +args.d+ "/FOvortex/" + str(int(time)) + "/"
		else:
			fPath = "./" +args.d+ "/FOvortex/" + str(time) + "/"

		for i in range(len(files)):
			fname = fPath + files[i]
			[cu,cr,ccf] = getSwirls(fname)
			qI = np.argmin(abs(qit - time))
			umag = uy[qI] + ccf
			Q = 3.14159*r2in*umag
			cs = 3.14159*cu*cr*cr/Q

			[miv,temp] = update(umag,miv,0,nOverall)
			[qiv,temp] = update(Q,qiv,0,nOverall)
			[mur,mur_std] = update(cu,mur,mur_std,nOverall)
			[mrr,mrr_std] = update(cr,mrr,mrr_std,nOverall)
			[msr,msr_std] = update(cs,msr,msr_std,nOverall)
			[tur,tur_std] = update(cu,tur,tur_std,nTime)
			[trr,trr_std] = update(cr,trr,trr_std,nTime)
			[tsr,tsr_std] = update(cs,tsr,tsr_std,nTime)
			nOverall += 1
			nTime += 1

		print("T = " + str(time) + " SR: " + str(round(tsr,2)) + " +/- " + str(round(tsr_std,4)) + "    Overall SR: " + str(round(msr,2)) + " +/- " + str(round(msr_std,4)))
		# sys.stdout.write('\r')
		# sys.stdout.write("[%-20s] %d%%" % ('='*i, 5*i))
		# sys.stdout.flush()


	print()
	print("Mean Inlet Velocity:  " + str(round(miv,4)) + " m/s")
	print("Mean Inlet Mass Flow: " + str(round(qiv,4)) + " m3/s")
	print("Mean Utheta Max:      " + str(round(mur,4)) + " m/s")
	print("Mean rUtheta Max:     " + str(round(mrr,4)) + " m")
	print("Mean Swirl Ratio:     " + str(round(msr,4)))
	print("Swirl Ratio StdDev:   " + str(round(msr_std,6)))
	print()

	return [mur,mrr]



if __name__ == "__main__":

	if checkArgs():

		if args.resid:
			residualPlot()

		if args.outflow:
			turn()

		if args.inlet:
			inletPlot()

		if args.profile or args.mat:
			vp()

		if args.calcSwirlRatio:
			csr()

		# [utmax, rutmax] = vortex()
		# uin = inlet(startIt,plot)
		# turn(startIt,plot)

		if not args.np:
			plt.show()