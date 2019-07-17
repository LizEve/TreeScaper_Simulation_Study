#!/bin/bash
from __future__ import division
import glob
import os
import dendropy as dp
import ast
import zipfile
import csv
import shutil
from collections import OrderedDict
import numpy
from sklearn.metrics.cluster import *
# We might try to divide by 0. its okay. 
numpy.seterr(divide='ignore', invalid='ignore')
# export PYTHONPATH=${PYTHONPATH}:'/usr/local/lib/python2.7/dist-packages'

### Unzipping scripts
def grabZip(zipFile,outFolder,suffix):
	zipZ = zipfile.ZipFile(zipFile, 'r')
	for file in zipZ.namelist():
		if file.endswith(suffix):
			zipZ.extract(file, outFolder)
			
def grabFiles(zipDir,outFolder,globber):
	os.chdir(zipDir)
	print(zipDir)
	print(globber)
	for f in glob.glob(globber):
		print(f)
		print("Unzip: "+str(f))
		#grabZip(f,outFolder,'Auto.out')
		#grabZip(f,outFolder,'distance.out')
		grabZip(f,outFolder,'Manu.out')
		grabZip(f,outFolder,'CovCommunities.out')
		#grabZip(f,outFolder,'Rec_Plateaus.out')
		grabZip(f,outFolder,'cloud.nex')
		
def grabNmove(zipDir,outDir,globber):
	# Define some paths
	globDir=os.path.join(outDir,globber.strip("results.zip").strip("_"))
	#print(outDir,globDir)
	# Get distance.out, Manu.out, CovCommunities.out, Plateaus.out, cloud.nex
	grabFiles(zipDir,outDir,globber)
	# Move files into main output folder
	os.chdir(globDir)
	os.system("mv * %s" % outDir)
	return outDir

# Make list of where to find starting trees in nexus file
def getStartingTreeNumList(nTrees,nCom):
	# make list to iterate through tree file and get starting trees
	nTreeList=[0]
	for x in range(0,nCom-2):
		n = nTreeList[x]+nTrees[x]
		nTreeList.append(n)
	return(nTreeList)


# Get simulated classes
# Returns: list -  list classes, each class a list of tree numbers, a dictionary of each tree number with associated class ID number. 
def getSimTrees(nCom,nTrees):
	# Create list of items(trees) in each simulated class
	# Create a dictionary with item:classID
	d = 0
	simComList=[]
	classIDdict={}
	for t in range(nCom-1):
		# Create numbered list of items for each sim class
		simCom = range(d,nTrees[t]+d)
		# Add itemIDs to list
		simComList.append(simCom)
		# Create dictionary with itemID and classID 
		for item in simCom:
			classIDdict[item]=t
		# Iterate through itemIDs as appropriate. 
		d = d + nTrees[t]
	return [simComList,classIDdict]

# Get found emp clusters
# Returns: list - number of found communities, list communities, each community a list of tree numbers, a dictionary of each tree number with associated cluster ID number. 
def getEmpTrees(affFile):
	# Read in analysis file
	file=open(affFile,'r').read()
	# Count number of clusters/communities found
	nEmpComA=file.count('Com')
	empComList=[]
	clusterIDdict={}
	# Read in list of items(trees) from each emp cluster
	for c in range(0,nEmpComA):
		for line in open(affFile,'r'):	
			# Grab string/line for each cluster
			if line.startswith("Community "+str(c+1)+" includes nodes: "):
					theline=line.replace("Community "+str(c+1)+" includes nodes: ","")
					# Change string to list
					theline= ast.literal_eval(theline)
					# Add itemIDs to list
					empComList.append(theline)
					# Create dictionary with itemID and clusterID 
					for item in theline:
						clusterIDdict[item]=c			
	return [nEmpComA,empComList,clusterIDdict]

# Set up class and cluster IDs in order
# Takes: two dictionaries. key=item value=class/clusterID. keys must be iterable
# Returns: two flat cluster/class ID lists 
def organizeIDLists(classIDdict,clusterIDdict):
	# Create two lists of class/cluster IDs that are ordered by shared itemIDs
	classIDList=[]
	clusterIDList=[]
	simnomatch=[]
	d=0
	f=0
	# different approach for type of key in dicts
	if isinstance(classIDdict.iterkeys().next(), int):
		# Sort dicts by key/itemID
		for sim in sorted(classIDdict.iterkeys()):
			for emp in sorted(clusterIDdict.iterkeys()):
				# If itemIDs match, add values (class/clusterID) to list
				if sim == emp:
					classIDList.append(classIDdict[sim])
					clusterIDList.append(clusterIDdict[emp])
	elif isinstance(classIDdict.iterkeys().next(), str):
		# Sort dicts by key/itemID
		for s in classIDdict.iterkeys():
			sim=set(ast.literal_eval(str(s)))
			d+=1
			for e in clusterIDdict.iterkeys():
				f+=1
				emp=set(ast.literal_eval(str(e)))
				# If itemIDs match, add values (class/clusterID) to list
				if sim == emp:
					classIDList.append(classIDdict[s])
					clusterIDList.append(clusterIDdict[e])
				else:
					simnomatch.append(s)
	else:
		print("ID dictionary keys are not a string or int, unsure how to proceed")
	# Check that class and cluster lists have the same number of items
	if len(classIDList)!=len(clusterIDList):
		print("unequal number of items")
	#print("check iterating")
	#print(d,f)
	return [classIDList,clusterIDList]

# NMI and other stats for class/cluster accuracy
# Returns a list of stats values
def accuracyCalcs(classIDList,clusterIDList):
	U=classIDList
	V=clusterIDList
	homo=homogeneity_score(U,V)
	comp=completeness_score(U,V)
	NMIsum=v_measure_score(U,V)
	NMIsqrt=normalized_mutual_info_score(U,V)
	aMI=adjusted_mutual_info_score(U,V)
	aRandInd=adjusted_rand_score(U,V)
	'''
	print("homogeneity: ",homo)
	print("completeness: ",comp)
	print("NMIsum: ",NMIsum)
	print("NMIsqrt: ",NMIsqrt)
	print("aMI: ",aMI)
	print("aRandInd: ",aRandInd)
	'''
	return[homo,comp,NMIsum,NMIsqrt,aMI,aRandInd]

# Get stats for Affinity
def runStatsA(affFile,nTrees,nCom,maxCom):
	# List of items(trees) in each simulated class
	# Dictionary of items with item:classID
	gst=getSimTrees(nCom,nTrees)
	#simComList = gst[0]
	classIDdict = gst[1]
	# Check if file exists
	if not os.path.isfile(affFile):
		# If it doesnt exist, just add lists of NA to fill spaces.
		classItemsA=nTrees
		totalClassA=len(classIDdict)
		affStats=[classItemsA,totalClassA]+["NA"]*9
	else:
		# List of items(trees) in each found emperical cluster
		# Dictionary with item:clusterID
		get=getEmpTrees(affFile)
		nClustA = get[0]
		empComList = get[1]
		clusterIDdict = get[2]

		# Create two lists of class/cluster IDs that are ordered by shared itemIDs
		oil=organizeIDLists(classIDdict,clusterIDdict)
		classIDList=oil[0]
		clusterIDList=oil[1]

		# Get information for output
		classItemsA=nTrees
		totalClassA=len(classIDList)
		clustItemsA=([len(x) for x in empComList])
		totalClustA=len(clusterIDList)
		acCalcsList=accuracyCalcs(classIDList,clusterIDList)
		homoA,compA,NMIsumA,NMIsqrtA,aMIA,aRandIndA = acCalcsList

		# Add everything to affstats
		affStats=[classItemsA,totalClassA,nClustA,clustItemsA,totalClustA,homoA,compA,NMIsumA,NMIsqrtA,aMIA,aRandIndA]
	print(affStats)
	return(affStats)

def getBipartTree(focalTree):
	bipartList = []
	for i in focalTree.internal_nodes():
		# Defines a list that initially contains all the leaf nodes from focal tree
		fullTaxonSet = focalTree.leaf_nodes()
		# Iterates over all internal nodes that are not the root
		if i is not focalTree.seed_node:
			# Instantiates string (conTree) to hold the constraint tree string
			conTree = "["
			# Iterates over leaf nodes that are descendants of the current internal node
			for j in i.leaf_nodes():
				# Appropriately adds the taxon name to the constraint tree string
				if j is i.leaf_nodes()[0]:
					conTree = conTree + str(j.taxon)
				else:
					conTree = conTree + "," + str(j.taxon)	
			# Closes out the part of the constraint for taxa descended from the focal node
			conTree = conTree + "]"
			# change from string to list
			conTree= ast.literal_eval(conTree)
			conTree=set(conTree)
			bipartList.append(conTree)
			#print conTree
	return bipartList

def getBipartComList(bipartFile,nCom):
	## Get community bipartitions
	# Start ordered dictionary 
	startStopDict = OrderedDict()
	# remember that nCom is desired coms + 1
	# create dictionary with stop/start lines for identifying communities in output file
	for c in range(1,nCom+1):
		start = 'Com '+str(c)+":"
		stop = 'Com '+str(c+1)+":"
		startStopDict[start]=stop
	# replace last value with a blank to match end of file
	last=startStopDict.keys()[-1]
	startStopDict[last]=''
	# Get dictionary of biparitions and their frequency in the community
	bipartComList = []
	for c in range(nCom):
		start=startStopDict.keys()[c]
		stop=startStopDict[start]
		#print(start,stop)
		comDict={}
		comList=[]
		with open(bipartFile,'r') as f:
			# Breaks the file into only the section we want. Found this online somehwere. 
			for line in f:
				if line.strip() == start:
					break
			for line in f:
				if line.strip() == stop:
					break
				# only print lines with lists
				if line.startswith("["):
					line = line.strip()
					# grab bipartition without trailing number
					bipart=line.split("]")[0]+"]"
					num=int((line.split("]")[1]).strip())
					# change from string to list
					bipart= ast.literal_eval(bipart)
					#print(bipart)
					#print(type(bipart))
					# Add list of bipartitions to list for community
					comDict[num] = set(bipart)
					comList.append(set(bipart))
		bipartComList.append(comList)
	return(bipartComList)

# Frequency of item in list
def freq(item,mylist,total):
	#print(mylist.count(item))
	frequency = mylist.count(item)/total
	return frequency

# List of unique items in list
def unqList(myList):
	# return unique items in a list of sets
	uniqueList = []
	for b in myList:
		# Use set for comparisons of bpart lists
		if set(b) not in uniqueList:
			uniqueList.append(b)
	return uniqueList

# Information on all trees/biparts in file
# Returns: nTreeList - list of where class breaks are in nexus file,
# allTrees - list of all trees in community
# fullBipList - list of all bipartitions in file. 
def getAllBiparts(nTrees,nCom,treeFile):
	# Get list of where starting trees are in file
	nTreeList=getStartingTreeNumList(nTrees,nCom)
	# Read in all tree strings from in file
	taxa = dp.TaxonNamespace()
	allTrees = dp.TreeList.get(path=treeFile,schema="nexus",rooting="default-rooted",taxon_namespace=taxa,preserve_underscores=True)
	#print(allTrees[50].as_string("newick"))
	# Get a list of all biparts in file
	fullBipList=[]
	# Get list of all bipartitions in file
	for t in allTrees:
		fullBipList.extend(getBipartTree(t))
	return [nTreeList,allTrees,fullBipList]

# List of bps that occur in 5-95% of trees in file
# Return - list of unique bipartitions
def filter595Biparts(totalTrees,allClassBps,unqBpList):
	# Get only biparts in >5% and <95% of trees
	unq595BpList=[]
	commonrarebiparts=[]
	for i in unqBpList:
		# Get frequency of bipartition in file
		frq = freq(i,allClassBps,totalTrees)
		#print(round(frq,3))
		if frq > 0.05 and frq < 0.95:
			unq595BpList.append(i)
		else:
			commonrarebiparts.append(i)
	return unq595BpList

# List of classes, each class a list of bipartitions, each bp a list of taxa
# Returns: List of classes with all bps found, and list of classes with bps that are unique and in 5-95% of trees
def getSimBiparts(nCom,nTreeList,allTrees,unq595BpList):
	# Get list of trees for each sim class, filtering for only biparts in 5-95% of trees in file, list of lists
	treeClassList=[]
	for t in range(0, nCom-1):
		if t != (nCom-2):
			clust = allTrees[nTreeList[t]:nTreeList[t+1]]
			treeClassList.append(clust)
		# special case for the last community go to end of file
		elif t == (nCom-2):
			clust = allTrees[nTreeList[t]:]
			treeClassList.append(clust)
	# Get list of bipartitions in each class
	simClassList=[]
	simClassAllList=[]
	# Iterate through classes of sim trees
	for c in treeClassList:
		classList=[]
		unqClassList=[]
		finalClassList=[]
		# For each tree in a class, get list of bipatitions
		for t in c:
			b=getBipartTree(t)
			classList.extend(b)
		unqClassList=unqList(classList)
		for b in unqClassList:
			if b in unq595BpList:
				finalClassList.append(b)
		simClassList.append(finalClassList)
		simClassAllList.append(classList)
	# returns list of all bps, and list of unique595 bps
	return [simClassAllList,simClassList]

# Returns list of classes with duplicates assigned to only one class
# If a bp found in both classes, the class where it is in higher frequency gets it. 
def filterDuplicates(simClassList,simClassAllList):
	# Get biparts shared between classes
	dupes=[]
	# iterate through each class
	for coms in range(len(simClassList)):
		# compare each bipartition in each class
		for bipart in simClassList[coms]:
			# to all other classes, but not self
			for com in range(len(simClassList)):
				if com != coms:
					if bipart in simClassList[com]:
						dupes.append(bipart)
	dupes=unqList(dupes)
	noDupSimClassList=[]
	for k in simClassList:
		noDupClassList=[]
		for bp in k:
			if bp not in dupes:
				noDupClassList.append(bp)
			elif bp in dupes:
				fracList=[]
				for ca in simClassAllList:
					x = ca.count(bp)
					frac = round((x/len(ca)),3)
					fracList.append(frac)
				majorFrac=max(fracList)
				majorClass=fracList.index(majorFrac)
				classNum=simClassList.index(k)
				#print(fracList)
				#print(majorFrac,majorClass)
				#print("test: "+str(classNum)+", "+str(majorClass))
				if majorClass == classNum:
					#print("add to: "+str(classNum)+", "+str(majorClass))
					noDupClassList.append(bp)
		noDupSimClassList.append(noDupClassList)
	return noDupSimClassList

# Turns list of lists into dictionary of items and IDs
# Returns: dictionary with keys as strings
def getIDdict(cList):
	cIDdict={}
	# Iterate through each cluster/class number
	for c in range(len(cList)):	
		# for each item in that cluster/class
		for i in cList[c]:
			# turn item into a list, then string
			x = str(list(i))
			# assign ID number to each item
			cIDdict[x]=c
	return(cIDdict)

# Get found emp clusters
# Returns: number of found communities, list - list communities, each community a list of bipartions
def getEmpBiparts(bipartFile):
	# Get number of emp communities found 
	file=open(bipartFile,'r').read()
	foundComs=file.count('Com')
	# Get total items in inferred clusters
	bipartComList=getBipartComList(bipartFile,foundComs)
	#print("Found clusters: "+str([foundComs,len(bipartComList)]))
	zz = []
	for xx in bipartComList:
		zz.append(len(xx))
	#print(zz,len(bipartComList[0])+len(bipartComList[1]))
	return [foundComs,bipartComList]

def runStatsB(treeFile, bipartFile, nTrees, nTips,nCom,maxCom):
	fileName=treeFile.strip(".nex")
	totalTrees=numpy.sum(nTrees)
	print(fileName)
	# Get information from simulated tree file
	gab=getAllBiparts(nTrees,nCom,treeFile)
	# List of class breakpoints
	nTreeList=gab[0]
	# List of all trees
	allTrees=gab[1]
	# List of all items(bipartitions) in file
	allClassBps=gab[2]
	# List of all unique bipatitions
	unqBpList = unqList(allClassBps)
	# List of biparts in >5% and <95% of trees
	unq595BpList=filter595Biparts(totalTrees,allClassBps,unqBpList)
	
	# List of items(bipartions) in each simulated class
	# filtering for only biparts in 5-95% of trees in file
	gsb=getSimBiparts(nCom,nTreeList,allTrees,unq595BpList)
	simClassAllList=gsb[0]
	simClassList=gsb[1]
	# Filter bps found in both classes, assign bp to class with higher occurance rate, turn into dictionary
	noDupSimClassList=filterDuplicates(simClassList,simClassAllList)
	classIDdict=getIDdict(noDupSimClassList)
	# if analysis didnt reach manu step. fill out with NA data
	if not os.path.isfile(bipartFile):
		classItemsC=([len(x) for x in noDupSimClassList])
		totalClassC=len(classIDdict)
		covStats=[classItemsC,totalClassC]+["NA"]*9
	else:
		# Get list of cluster biparts, list of lists of sets, into dict
		geb=getEmpBiparts(bipartFile)
		nClustC = geb[0]
		bipartComList=geb[1]
		clusterIDdict=getIDdict(bipartComList)
		#print("len on bipart lists")
		#print(len(noDupSimClassList[0]),len(noDupSimClassList[1]),len(bipartComList[0]),len(bipartComList[1]))
		#print("len of dicts")
		#print(len(classIDdict),len(clusterIDdict))
		# Create two lists of class/cluster IDs that are ordered by shared itemIDs
		oil=organizeIDLists(classIDdict,clusterIDdict)
		classIDList=oil[0]
		clusterIDList=oil[1]
		#print("ID list info: ")
		#print(classIDList)
		#print(clusterIDList)
		#print(len(classIDList))
		#print(len(clusterIDList))
		# Get information for output
		classItemsC=([len(x) for x in noDupSimClassList])
		totalClassC=len(classIDList)
		clustItemsC=([len(x) for x in bipartComList])
		totalClustC=len(clusterIDList)
		acCalcsList=accuracyCalcs(classIDList,clusterIDList)
		homoC,compC,NMIsumC,NMIsqrtC,aMIC,aRandIndC = acCalcsList
		# Add everything to covstats
		covStats=[classItemsC,totalClassC,nClustC,clustItemsC,totalClustC,homoC,compC,NMIsumC,NMIsqrtC,aMIC,aRandIndC]
	print(covStats)
	return covStats

# Get number of trees in each community
# Returns: list 
def parseNTrees(nTree,nCom):
	# parse trees for asymetrical communities, for multicom the first cluster in the nexus file is always the smallest.
	if "v" not in nTree:
		nTrees = []
		for n in range(0,nCom-1):
			nTrees.append(int(nTree))
		#print(nTrees)
	elif "v" in nTree:
		# for asym com, use #v# and total(which is not the actual total for multiple communities) to get number of trees per cluster
		com1v = float(nTree.split("v")[0])
		com2v = float(nTree.split("v")[1])
		total = 1000
		com1nTrees= int(total*round((com1v/(com1v+com2v)),3))
		com2nTrees= int(total*round((com2v/(com1v+com2v)),3))
		# if there are multiple coms, add multiple big coms
		if nCom > 3:
			nTrees = [com1nTrees]
			for n in range(1,nCom-1):
				nTrees.append(com2nTrees)
		elif nCom == 3:
			smallFirst=[2,100,500,1000]
			bigFirst=[5,10,20,50]
			if int(com2v) in smallFirst:
				nTrees=[com1nTrees,com2nTrees]
			elif int(com2v) in bigFirst:
				nTrees=[com2nTrees,com1nTrees]
			else:
				print("I'm sorry dave, I'm afraid I cant do that")
	return nTrees

def singleFileStats(treeFile,nCom,maxCom):
	nCom = nCom + 1
	treeSetTrunc=treeFile.strip(".nex")
	# Get important stats from file name
	nTips = float(treeSetTrunc.split("_")[0].strip("tip"))
	nTree = treeSetTrunc.split("_")[1].strip("trees")
	nTrees = parseNTrees(nTree,nCom)
	# Get files
	bipartFile=treeSetTrunc+"_CPM_CovCommunities.out"
	affFile=treeSetTrunc+"_CPM_Rec_AffManu.out"
	# Get biparition community stats
	statsB = runStatsB(treeFile,bipartFile,nTrees,nTips,nCom,maxCom)
	# Get affinity community stats
	statsA = runStatsA(affFile,nTrees,nCom,maxCom)
	# for testing only
	#naList = ["NA"]*16
	#statsB = naList
	#statsA = naList
	statsAll = [treeSetTrunc,nTips,nCom-1]
	statsAll.extend(statsA)
	statsAll.extend(statsB)
	print(statsAll)
	return(statsAll)

def inOutResults(outDir,outDataFile,nCom,maxCom):
	headerLine = ["FileName","nTips","nSimCom","classItemsA","totalClassA","nClustA","clustItemsA","totalClustA","homoA","compA","NMIsumA","NMIsqrtA","aMIA","aRandIndA","classItemsC","totalClassC","nClustC","clustItemsC","totalClustC","homoC","compC","NMIsumC","NMIsqrtC","aMIC","aRandIndC"]
	os.chdir(outDir)
	#Open and prep dataSheet
	print("Gathering data for: "+str(outDataFile)+" from "+str(outDir))
	with open(outDataFile, 'a') as fp:
		wr = csv.writer(fp, dialect='excel')
		wr.writerow(headerLine)
		# Iterate through all files and folders
		#for f in glob.glob("*.nex"):
		for f in glob.glob("*cloud"):
			os.chdir(f)
			#print("##########################")
			#print(f)
			# grab data
			data = singleFileStats(f+".nex",nCom,maxCom)
			#print(len(headerLine))
			#print(len(data))
			# write data
			wr.writerow(data)
			os.chdir(outDir)

def iterateMain(myList,nCom,zipDir,maxCom):
	for zipFile in myList:
		folderName=zipFile.strip("results.zip").strip("_")
		globDir=os.path.join(zipDir,folderName)
		outDataFile=os.path.join(zipDir,folderName+'.csv')
		os.chdir(zipDir)
		grabFiles(zipDir,zipDir,zipFile)
		os.chdir(globDir)
		inOutResults(globDir,outDataFile,nCom,maxCom)


### MAIN
def main():
	# Constants
	# Maximum number of communities specified
	maxCom = 10
	# Number of communities simulated
	nCom = 2
	# Folder where zip file of results is
	zipDir='/media/ExtraDrive4/gmount/TreeScaperResults'
	# Name of the zip file
	zipFile='numTips100_results.zip'
	# Name of the folder when uncipped
	globDir=os.path.join(zipDir,'numTips100')
	# File to put results in
	outDataFile=os.path.join(zipDir,'numTips100.csv')
	os.chdir(zipDir)
	grabFiles(zipDir,zipDir,zipFile)
	os.chdir(globDir)
	inOutResults(globDir,outDataFile,nCom,maxCom)


if __name__=='__main__':
	main()


