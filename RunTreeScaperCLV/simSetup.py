#! /usr/bin/env python

import os
import glob
import shutil 
import itertools
import re 


# Make sure you are in directory that you are calling script from 
mainDir = os.getcwd()

for t in glob.glob('*.nex'):

	# Get base name of file
	treeSet=str(t)
	treeSetIndex = treeSet.find(".nex")
	fName = treeSet[:treeSetIndex]
	
	# Create path 
	dirPath = os.path.join(mainDir,fName)
	print(dirPath)

	# Make directory 
	if not os.path.exists(dirPath):
		#print("NEW DIR")
		os.mkdir(dirPath)	

	# Copy nexus to new folder
	os.system("cp %s %s" % (t,dirPath))
