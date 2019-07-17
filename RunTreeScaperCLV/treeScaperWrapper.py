#!/usr/bin/env python
# -*- coding: utf-8 -*-


#### Output files
# _Cov/AffPlateaus.out -> Matrix of each plateau found, length, lower and upper bounds, mean lambda, and number of plateaus. Ordered by largest to smallest plateau size
# _Cov/AffAuto.out -> log file of auto plateau finder run

import re
import os
import sys
import glob
import numpy as np
import time
import dendropy
from shutil import move, copyfile
from tempfile import mkstemp
import argparse 
import operator
import csv


__author__ = 'Genevieve G. Mount, Jeremy Ash'
__version__ = 'July 2019'

def get_args():
	'''This function parses and return arguments passed in'''
	parser = argparse.ArgumentParser("Wrapper for CLVTreeScaper")
	# Add arguments
	parser.add_argument("-c", "--clvPath", help="Path to your command line treescaper, path/to/CLVTreeScaper", required=True)
	parser.add_argument("-i", "--inputNexus", help="Nexus file with translate block. This file will be modified. Trees will be numbered [1]-[#trees]", required=True)
	parser.add_argument("-n", "--network", help="Type of network used to identify communities, can be Covariance/Affinity",required=True)
	parser.add_argument("-f", "--findPlateau",help="Find plateau automatically, manually, or automatically and then run manually, Auto/Manual/Both",required=True)
	parser.add_argument("-m", "--model", help="Model for community detection, can be CNM/CPM(default)/ERNM/NNM",required=False, default='CPM')
	parser.add_argument("-w", "--weighted", help="0 = unweighted(default), no branch lengths, 1 = weighted, branch lengths", required=False, default=0, type=int)
	parser.add_argument("-r", "--rooted",help="0=unrooted, 1=rooted(default)",required=False, default=1, type=int)
	parser.add_argument("-lf", "--lowFreq",help="For covariance networks, nodes with frequencies below this value are ignored. A number between 0 and 1 (default=0.05)",required=False, default=0.05, type=float)
	parser.add_argument("-hf", "--highFreq",help="For covariance networks, nodes with frequencies above this value are ignored. A number between 0 and 1 (default=0.95)",required=False, default=0.95, type=float)
	parser.add_argument("-dm", "--distanceMatrix",help="For affinity network. Indicates the distance metric. Options- ‘URF’: Unweighted Robinson-Foulds (default)/ ‘RF’: Weighted Robinson-Foulds/ ‘Mat’: Matching distance/ ‘SPR’: Subtree-Prune-Regraft.",required=False, default="URF")
	parser.add_argument("-at", "--affinityTransformation",help="For affinity network. Indicates the type of distand to affinity transformation. Options- ‘Rec’: Reciprocal(default). ‘Exp’: Exponential",required=False, default='Rec')
	parser.add_argument("-p", "--plateau",help="Input absolute lambda value from desired plateau, then indicate pos/neg. ex: '0.3,neg' would mean -0.3", required=False, default="5000000,placeholder", type=str)
	parser.add_argument("-mc", "--maxCommunityNumber",help="Maximum number of communities reasonable. This avoids oversplitting. Default=100",required=False, default=100, type=int)
	parser.add_argument("-s", "--simulationStudy", help="For simulation, study, see below for more info.", required=False, default="X,1:1,X")
	#### Simulation extras, need placeholders(X) for options not being used.
	# -s "gnu" - jumps in and out of folder based on input nexus name. make sure proper CLV path is specified. Often will be -c "../CLVTreeScaper"
	# -s "X,treeRatio" - ratio of size of clusters. MUST input in order that clouds are in file. 
	# -s "X,1:1,fff" - FindFromFile, find largest plateau from already run auto plateau run. input info and filename must match exactly

	# Array for all arguments passed to script
	args = parser.parse_args()
	# Assign args to variables
	clvPath = args.clvPath
	inNexus = args.inputNexus
	findPlateau = args.findPlateau
	model = args.model
	network = args.network
	rooted = args.rooted
	weighted = args.weighted
	plateau = float(args.plateau.split(",")[0])
	platSign = args.plateau.split(",")[1].strip()
	maxCom = args.maxCommunityNumber
	# Covariance
	hf = args.highFreq
	lf = args.lowFreq
	# Affinity
	dm = args.distanceMatrix
	at = args.affinityTransformation
	# Set lambda negative to zero. lambda pos is the abs value of a lambda within the chosen plateau
	lamNeg = 0
	sim = args.simulationStudy.split(",")
	# Grabbing naming info
	treeSet=str(inNexus)
	treeSetIndex = treeSet.find(".nex")
	treeSetTrunc = treeSet[:treeSetIndex]
	if network == 'Affinity':
		outLog = treeSetTrunc+"_"+model+"_"+at
	elif network == 'Covariance':
		outLog = treeSetTrunc+"_"+model
	# Return all variable values
	return clvPath, model, findPlateau, network, rooted, weighted, plateau, platSign, maxCom, hf, lf, dm, at, lamNeg, sim, treeSet, treeSetTrunc, outLog

def nRF(RF,tips):
	maxRF = 2*(float(tips)-2)
	norm = float(RF)/maxRF
	return round(norm, 5)

def make_list(pre_ls,convert):
	# Turns tab delimited line into a list. 
	pre_ls = pre_ls.split("\t")
	ls = [i for i in pre_ls if i != "\n"]
	if convert == "int":
		ls = [int(i) for i in ls]
	if convert == "float":
		ls = [float(i) for i in ls]
	#ls = ls[1:]
	#print("list - "+str(ls))
	return ls

def list_to_file(header, majorList, outFileName):
	# Takes list, header, and outfile
	# Writes out a file with header and lists. 
	with open(outFileName, "w") as outFile:
		outFile.write("%s\n" % header)
		wr = csv.writer(outFile)
		wr.writerows(majorList)

def reg_ex_match(file, pattern):
	# Returns the first match of a reg ex search
	file.seek(0)
	for line in file:
		m = pattern.match(line)
		if m:
			return m.group(1)

def edit_treeset(treeFileEditPath):
	# Adds comment blocks that number each tree with the indices used by TreeScaper. 
	treeFileEdit = open(treeFileEditPath, 'r')

	lineNum = 1
	# make a temp file
	fh, absPath = mkstemp()
	tempFile = open(absPath,'w')
	for line in treeFileEdit:
		if line.find('[&U]') != -1: # Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
			if line.find('['+str(lineNum)+']') != -1:
				tempFile.write(line)
			else:
				tempFile.write(line.replace('=','['+str(lineNum)+']='))
			lineNum += 1
		elif line.find('[&R]') != -1: # Need to have this in every line with a tree! This will be [&U] for unrooted trees and [&R] for rooted
			if line.find('['+str(lineNum)+']') != -1:
				tempFile.write(line)
			else:
				tempFile.write(line.replace('=','['+str(lineNum)+']='))
			lineNum += 1
		else:
			tempFile.write(line)
	# close temp file
	tempFile.close()
	os.close(fh)
	treeFileEdit.close()
	# Remove original file
	os.remove(treeFileEditPath)
	# Move new file
	move(absPath, treeFileEditPath)
	return lineNum

def time_calc(startTime, endTime):
	totalTime = str(round(endTime - startTime, 5))
	return totalTime

def time_calc_both(startTimeAuto,endTimeAuto,startTimeFindL,endTimeFindL,startTimeManu,endTimeManu,startTimeParseOut, endTimeParseOut):
		autoTime = time_calc(startTimeAuto,endTimeAuto)
		findLTime = time_calc(startTimeFindL,endTimeFindL)
		manuTime = time_calc(startTimeManu,endTimeManu)
		parseOutTime = time_calc(startTimeParseOut, endTimeParseOut)
		return [autoTime, findLTime, manuTime, parseOutTime]

def get_complicated_filename(treeSetTrunc,rooted,weighted,dm):
	if dm == 'Mat':
		dmM = 'Matching'
	else:
		dmM = dm
	if rooted == 1:
		if weighted == 0:
			filename=treeSetTrunc+"_rooted_unweighted_"+dmM+"-distance.out"
		elif weighted == 1:
			filename=treeSetTrunc+"_rooted_weighted_"+dmM+"-distance.out"
	elif rooted == 0:
		if weighted == 0:
			filename=treeSetTrunc+"_unrooted_unweighted_"+dmM+"-distance.out"
		elif weighted == 1:
			filename=treeSetTrunc+"_unrooted_weighted_"+dmM+"-distance.out"
	if dm != 'URF':
		print("there may(or may not) be problems reading in distance.out file.\nEdit function get_complicated_filename if needed")
	return filename

def output_crap(treeSetTrunc, weighted, rooted, network, model, dm, at, hf, lf, lamNeg, maxCom, platSign, plateauInfo, timing, sim, outLog):
	# Grab a bunch of info to output for simulation study. Implemented in "Both" plateau finder option. 
	tips = treeSetTrunc.split("_")[0].strip("tip")
	trees = treeSetTrunc.split("_")[1].strip("trees")
	cloudRatio = sim[1]
	cloudDensity = treeSetTrunc.split("_")[2]
	startRF = treeSetTrunc.split("_")[3].strip("start")
	edgeRF = float(startRF)-(float(cloudDensity)*2)
	if plateauInfo[0] == "NA":
		plateauFound, numComsUsed, orderFound = "NA", "NA", "NA"
	else: 
		plateauFound = abs(plateauInfo[0])
		numComsUsed = plateauInfo[1]
		orderFound = plateauInfo[2]

	allComsFound = plateauInfo[3][:6]
	orderComsFound = plateauInfo[4][:6]
	lenComsFound = plateauInfo[5][:6]
	# Pull out top 4 plateaus and their info. Populate lists with NA if there are less than 4 communities found. 
	numComsList=["NA", "NA","NA", "NA", "NA", "NA"]
	orderComsList=["NA", "NA","NA", "NA", "NA", "NA"]
	lenComsList=["NA", "NA","NA", "NA", "NA", "NA"]
	count = 0
	for num in allComsFound:
		numComsList[count]=num
		count+=1
	count = 0
	for num in orderComsFound:
		orderComsList[count]=num
		count+=1
	count = 0
	for num in lenComsFound:
		lenComsList[count]=num
		count+=1
	#print(numComsList)
	#print(orderComsList)
	#print(lenComsList)
	autoTime = timing[0]
	findLTime = timing[1]
	manuTime = timing[2]
	parseOutTime = timing[3]
	totalTime = float(autoTime) + float(findLTime) + float(manuTime) + float(parseOutTime)
	if network == 'Covariance':
		dm,at = "NA","NA"
	if network == 'Affinity':
		hf,lf = "NA","NA"
	header = ("FileName,Tips,Trees,Asymmetry,StartDistance,EdgeDistance,Weighted,Rooted,Network,Model,DistanceMatrix,AffinityTransformation,HighFreq,LowFreq,LambdaNegative,MaximumCommunitiesAllowed,LambdaSign,LambdaUsed,NumberCommunities,OrderFound,oneNumCom,oneComOrder,oneComLen,twoNumCom,twoComOrder,twoComLen,threeNumCom,threeComOrder,threeComLen,fourNumCom,fourComOrder,fourComLen,fiveNumCom,fiveComOrder,fiveComLen,sixNumCom,sixComOrder,sixComLen,TimeTotal,TimeManual,TimeAuto,TimeFindLambda,TimeParseOutput\n")
	outFile = open( "%s_%s_log.out" % (outLog, network) , 'w' )
	outFile.write(header)
	outFile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (treeSetTrunc,tips,trees,cloudRatio,startRF,edgeRF,weighted,rooted,network,model,dm,at,hf,lf,lamNeg,maxCom,platSign,plateauFound,numComsUsed,orderFound,numComsList[0],orderComsList[0],lenComsList[0],numComsList[1],orderComsList[1],lenComsList[1],numComsList[2],orderComsList[2],lenComsList[2],numComsList[3],orderComsList[3],lenComsList[3],numComsList[4],orderComsList[4],lenComsList[4],numComsList[5],orderComsList[5],lenComsList[5],totalTime,autoTime,manuTime,findLTime,parseOutTime))
	print("lambda used")
	print(plateauFound)
	print("numcomchosen")
	print(numComsUsed)
	print("numComsList")
	print(numComsList)
	print("len coms found")
	print(lenComsList)

def auto_find_lambda(inOutFile):
	# Pulls all plateaus from auto output file. 
	# Returns list of lists sorted with largest plateau first.
	# Each list for each plateau are in the format:
	# [length, lowerBound, upperBound, mean labda value of plateau, number of communities]
	pattern = 'Updated plateau:'
	plateaus=[]
	with open(inOutFile , 'r' ) as file:
		for line in file:
			if pattern in line:
				# Store lines in a list to pull out more info. 
				plateaus.append(line)
	majorList = []
	for x in plateaus:
		# Pull out line with updated plateau
		lineParsed = x.split(":")
		# Do some python gymnastics to pull out just the numbers
		lowerBound = round(float(str(lineParsed[1].split(",")[0]).strip().strip("[")), 6)
		upperBound = round(float(str(lineParsed[1].split(",")[1]).strip().strip("]")), 6)
		length = round(float(str(lineParsed[2].split(",")[0]).strip()), 6)
		numCom = int(str(lineParsed[3]).strip())
		# Get lambda value from middle of plateau
		plat = [lowerBound, upperBound]
		manLambda = round(np.mean(plat), 4)
		# Or get middle of 0 and upper bound
		if manLambda == 0:
			manLambda = round(np.mean([0,plat[1]]), 4)
		# Create list of important values
		minorList = [length, lowerBound, upperBound, manLambda, numCom]
		majorList.append(minorList)

	# Sort list by size of plateau
	majorList.sort(key=operator.itemgetter(0), reverse=True)
	return majorList

def limit_num_communities(majorList, maxCom):
	# Returns first instance of lambda that corresponds to a plateau with less than maxCom communities. 
	# This first instance should be the largest plateau with num com under given maxCom
	# majorList = list for each plateau [plateaulength, upperbound, lowerbound, meanlambda, num coms]
	#print("[plateaulength, upperbound, lowerbound, meanlambda, num coms]")
	#print(majorList)
	numPlat = (len(majorList))
	orderComsFound=[]
	allComsFound=[]
	lenComsFound=[]
	comDict={}
	plateauFound="NA"
	numComsUsed="NA"
	orderFound="NA"
	# Gather information about all communities
	for plat in range(numPlat):
		allComsFound.append(majorList[plat][4])
		orderComsFound.append(plat+1)
		lenComsFound.append(majorList[plat][0])
	for plat in range(numPlat):
		numcoms = majorList[plat][4]
		if numcoms <= maxCom:
			# If there are less than maxCom communities
			# make dictionary: plat len: numcom, mean lambda, order,
			plateauSize=majorList[plat][0]
			meanLambda=majorList[plat][3]
			numComs=majorList[plat][4]
			orderFound=plat+1
			comDict[plateauSize]=[numComs,meanLambda,orderFound]
	if len(comDict) != 0:
		biggest=max(comDict)
		numComsUsed=comDict[biggest][0]
		plateauFound=comDict[biggest][1]
		orderFound=comDict[biggest][2]
	# return meanLambda for chosen plateau, number of communites, order in which this plateau was found, list of all communities found
	# plateauFound, numComsUsed, orderFound, allComsFound
	return[plateauFound, numComsUsed, orderFound,allComsFound,orderComsFound,lenComsFound]

def parse_cov_output(treeSetTrunc, model):
	# Takes manual covariance output files *_CovManu.out 
	# Pulls out community and reports which bipartitions are in each community
	# Reports number of trees that contain each biparition
	# Returns number of communities
	outLog = treeSetTrunc+"_"+model
	# Open out file from manual run
	comFile = open("%s_CovManu.out" %  outLog, 'r')
	pattern = re.compile('Number of communities: (\d+)')
	coms = int(reg_ex_match(comFile, pattern))
	comKey = open("%s_CovCommunities.out" %  outLog, 'w')
	# Makes a key to decode the bipartitions in each community
	for i in range(1, coms+1): 
		# Get nodes for each community
		pattern = re.compile('Community '+str(i)+' includes nodes: (.+)')
		comStr = reg_ex_match(comFile, pattern)
		comLS = comStr.split(",")
		print('Community '+str(i)+' includes nodes: '+str(comLS))
		comLS = filter(None, comLS)
		comKey.write("Com %s:\n" % i)
		# For each node
		for j in comLS:
			# Bipartitions count from 1 while nodes count from 0.
			j = int(j)+1
			pattern = re.compile("bipartition "+str(j)+" : ([0-1]+?), appear times: ([0-9]+?)$")
			# Some file copying, temporary
			fh, absPath = mkstemp()
			copyfile("%s_CovManu.out" %  outLog, absPath)
			comTempFile = open(absPath,'r')
			# Grab bipartitions
			for line in comTempFile:
				m = pattern.match(line)
				if m:
					# Grab binary code for bipartition
					bipart = m.group(1)
					print("Bipartion "+str(j)+" : "+str(bipart))
					# Grab number of times bipartition occurs
					freq = m.group(2)
					bipartLS = []
					# Create list from binary code
					for c in bipart:
						bipartLS.append(c)
					# Create a list of taxa included in each bipartition
					indices = [x+1 for x, y in enumerate(bipartLS) if y == '1']
					# Get taxon names
					for k in indices:
						pattern2 = re.compile("(.+) , "+str(k))
						taxon = reg_ex_match(comFile, pattern2)
						indices[indices.index(k)] = taxon
			# Close temp file
			comTempFile.close()
			os.close(fh)
			# Write taxon and frequency for bipartitions in community
			comKey.write("%s %s\n" % (indices,freq))
		print("\n")
		comKey.write("\n")
	comKey.close()
	return coms

def run_covariance_auto(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, hf, lf):
	outLog = treeSetTrunc+"_"+model
	os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Covariance -cm %s -lm auto -hf %s -lf %s > %s_CovAuto.out" % (clvPath, treeSet, weighted, rooted, model, hf, lf, outLog))

def run_covariance_man(clvPath, treeSet, treeSetTrunc, weighted, rooted, model,  hf, lf, lamNeg, plateauFound, platSign):
	outLog = treeSetTrunc+"_"+model
	if platSign == 'pos':
			os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Covariance -cm %s -lm manu -lp %s -ln %s -hf %s -lf %s > %s_CovManu.out" % (clvPath, treeSet, weighted, rooted, model, plateauFound, lamNeg, hf, lf, outLog))
	# Avoid passing negative number to CLV
	# Pass absolute value of lambda to lambda negative and 0 to lambda positive
	elif platSign == 'neg':
		plateauFound = abs(plateauFound)
		os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Covariance -cm %s -lm manu -lp %s -ln %s -hf %s -lf %s > %s_CovManu.out" % (clvPath, treeSet, weighted, rooted, model, lamNeg,plateauFound, hf, lf, outLog))
	else:
		print("Indicate 'pos' or 'neg' for the -p --Plateau flag. Ex: -p 0.05,'neg'")

def run_covariance_both(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, hf, lf, lamNeg, maxCom, outLog, sim):
	# zero out platSign as safety measure
	platSign=None
	if sim[2] == "X":
		print("Running auto")
		# Run automatic plateau finder
		startTimeAuto = time.time()
		run_covariance_auto(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, hf, lf)
		endTimeAuto = time.time()
	elif sim[2] == "fff":
		print("not running auto")
		startTimeAuto = 0
		endTimeAuto = 0
	else:
		print("Sim arguments are wrong. Quitting...")
		raise SystemExit
	# Find all plateaus
	startTimeFindL = time.time()
	majorList = auto_find_lambda("%s_CovAuto.out" % outLog)
	# Output plateau info to file
	header = "length, lowerBound, upperBound, meanLambda, numCom"
	outFileName = ("%s_CovPlateaus.out" % outLog)
	list_to_file(header, majorList, outFileName)
	# Limit number of communities identified
	# Save info: plateauFound, numComsUsed, orderFound, allComsUsed[list]
	plateauInfo = limit_num_communities(majorList, maxCom)
	endTimeFindL = time.time()
	# Pull out plateau values or print output and quit
	try:
		int(plateauInfo[0])
		plateauFound = plateauInfo[0]
	except ValueError:
		print("\n\nThere are no plateaus that have "+str(maxCom)+" or less communities. \nTry visualizing plateaus using GUI. \nOutputting information and then quitting...\n\n")
		timing = time_calc_both(startTimeAuto,endTimeAuto,startTimeFindL,endTimeFindL,1,1,1,1)
		output_crap(treeSetTrunc, weighted, rooted, network, model, dm, at, hf, lf, lamNeg, maxCom, "NA", plateauInfo, timing, sim, outLog)
		raise SystemExit
	# Assess sign of plateau lambda
	if plateauFound > 0:
		platSign='pos'
	elif plateauFound < 0:
		platSign='neg'
	else:
		print("\n\nSomething has gone very wrong.\nQuitting...\n\n")
		raise SystemExit
	# Run manual with found plateau
	startTimeManu = time.time()
	run_covariance_man(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, hf, lf, lamNeg, plateauFound, platSign)
	endTimeManu = time.time()
	# Parse the output from the manual run
	startTimeParseOut = time.time()
	parse_cov_output(treeSetTrunc, model)
	endTimeParseOut = time.time()
	# Record timing of all parts
	# returns a list of [autoTime, findLTime, manuTime, parseOutTime]
	timing = time_calc_both(startTimeAuto,endTimeAuto,startTimeFindL,endTimeFindL,startTimeManu,endTimeManu,startTimeParseOut, endTimeParseOut)
	# output run info to file
	output_crap(treeSetTrunc, weighted, rooted, network, model, dm, at, hf, lf, lamNeg, maxCom, platSign, plateauInfo, timing, sim, outLog)

def run_affinity_auto(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, outLog):
	os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Affinity -cm %s -lm auto -dm %s -am %s > %s_AffAuto.out" % (clvPath, treeSet, weighted, rooted, model, dm, at, outLog))

def run_affinity_man(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, lamNeg, plateauFound, platSign, outLog):
	if platSign == 'pos':
		os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Affinity -cm %s -lm manu -dm %s -am %s -lp %s -ln %s > %s_AffManu.out" % (	clvPath, treeSet, weighted, rooted, model, dm, at, plateauFound, lamNeg, outLog))
	# Avoid passing negative number to CLV
	# Pass absolute value of lambda to lambda negative and 0 to lambda positive
	elif platSign == 'neg':
		plateauFound = abs(plateauFound)
		os.system("%s -trees -f %s -ft Trees -w %s -r %s -o Community -t Affinity -cm %s -lm manu -dm %s -am %s -lp %s -ln %s > %s_AffManu.out" % (clvPath, treeSet, weighted, rooted, model, dm, at, lamNeg, plateauFound, outLog))
	else:
		print("Indicate 'pos' or 'neg' for the -p --Plateau flag. Ex: -p 0.05,'neg'")

def run_affinity_both(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, lamNeg, maxCom, outLog, sim):
		# zero out platSign as safety measure
		platSign=None
		if sim[2] == "X":
			print("Running auto")
			# Run automatic plateau finder
			startTimeAuto = time.time()
			run_affinity_auto(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, outLog)
			endTimeAuto = time.time()
		elif sim[2] == "fff":
			print("not running auto")
			startTimeAuto = 0
			endTimeAuto = 0
		else:
			print("Sim arguments are wrong. Quitting...")
			raise SystemExit

		# Print out all plateaus
		startTimeFindL = time.time()
		majorList = auto_find_lambda("%s_AffAuto.out" % outLog)

		# Output plateau info to file
		header = "length, lowerBound, upperBound, meanLambda, numCom"
		outFileName = ("%s_AffPlateaus.out" % outLog)
		list_to_file(header, majorList, outFileName)

		# Limit number of communities identified
		# Save info: plateauFound, numComsUsed, orderFound, allComsUsed[list]
		plateauInfo = limit_num_communities(majorList, maxCom)
		endTimeFindL = time.time()

		# Pull out plateau values or print output and quit
		try:
			int(plateauInfo[0])
			plateauFound = plateauInfo[0]
		except ValueError:
			print("\n\nThere are no plateaus that have "+str(maxCom)+" or less communities. \nPlease adjust this number and try again. \nProgram outputting information and then quitting...\n\n")
			timing = time_calc_both(startTimeAuto,endTimeAuto,startTimeFindL,endTimeFindL,1,1,1,1)
			output_crap(treeSetTrunc, weighted, rooted, network, model, dm, at, hf, lf, lamNeg, maxCom, "NA", plateauInfo, timing, sim, outLog)
			raise SystemExit
		if plateauFound > 0:
			platSign='pos'
		elif plateauFound < 0:
			platSign='neg'
		else:
			print("\n\nSomething has gone very wrong.\nQuitting...\n\n")
			raise SystemExit
		# Run manual with found plateau
		startTimeManu = time.time()
		run_affinity_man(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, lamNeg, plateauFound, platSign, outLog)
		endTimeManu = time.time()
		# Parse the output from the manual run
		startTimeParseOut = time.time()
		endTimeParseOut = time.time()
		# Record timing of all parts
		# returns a list of [autoTime, findLTime, manuTime, parseOutTime]
		timing = time_calc_both(startTimeAuto,endTimeAuto,startTimeFindL,endTimeFindL,startTimeManu,endTimeManu,startTimeParseOut, endTimeParseOut)
		# output run info to file
		output_crap(treeSetTrunc, weighted, rooted, network, model, dm, at, hf, lf, lamNeg, maxCom, platSign, plateauInfo, timing, sim, outLog)

def main():
	# Match return values from get_arguments()
	# and assign to their respective variables
	clvPath, model, findPlateau, network, rooted, weighted, plateau, platSign, maxCom, hf, lf, dm, at, lamNeg, sim, treeSet, treeSetTrunc, outLog = get_args()
	try:
		file = open(treeSet, 'r')
	except IOError:
		print('\n\nThere was an error opening the file.\nPlease check filename.\nQuitting...\n\n')
		raise SystemExit

	if network == 'Affinity':
		if findPlateau == 'Auto':
			run_affinity_auto(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, outLog)
			# Print out all plateaus
			majorList = auto_find_lambda("%s_AffAuto.out" % outLog)
			header = "length, lowerBound, upperBound, meanLambda, numCom"
			outFileName = ("%s_AffPlateaus.out" % outLog)
			list_to_file(header, majorList, outFileName)

		elif findPlateau == 'Manual':
			run_affinity_man(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, lamNeg, plateau, platSign, outLog)
		
		elif findPlateau == 'Both':
			run_affinity_both(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, dm, at, lamNeg, maxCom, outLog, sim)

		else:
			print("Please enter option for --findPlateau, 'Auto', 'Manual' or 'Both'")

	elif network == 'Covariance':
		if findPlateau == 'Auto':
			run_covariance_auto(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, hf, lf)
			# Print out all plateaus
			majorList = auto_find_lambda("%s_CovAuto.out" % outLog)
			header = "length, lowerBound, upperBound, meanLambda, numCom"
			outFileName = ("%s_CovPlateaus.out" % outLog)
			list_to_file(header, majorList, outFileName)

		elif findPlateau == 'Manual':
			run_covariance_man(clvPath, treeSet, treeSetTrunc, weighted, rooted, model,  hf, lf, lamNeg, plateau, platSign)
			parse_cov_output(treeSetTrunc, model)
		
		elif findPlateau == 'Both':
			run_covariance_both(clvPath, treeSet, treeSetTrunc, weighted, rooted, model, hf, lf, lamNeg, maxCom, outLog, sim)
		else:
			print("Please properly enter option for --findPlateau, 'Auto', 'Manual' or 'Both'")


if __name__=='__main__':
	clvPath, model, findPlateau, network, rooted, weighted, plateau, platSign, maxCom, hf, lf, dm, at, lamNeg, sim, treeSet, treeSetTrunc, outLog = get_args()
	# If running with gnu parallel, switch into each loci folder, perform action, and then come back out. 
	if sim[0] == "gnu":
		mainDir = os.getcwd()
		dirPath = os.path.join(mainDir,treeSetTrunc)
		os.chdir(dirPath)
		print("Changing to folder: "+str(dirPath))
		main()
		os.chdir(mainDir)
	else:
		main()