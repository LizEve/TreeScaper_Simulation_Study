#! /usr/bin/env python
import os
import glob
import subprocess

# File name example
# 10tip_1000trees_0.125_0.5start_10_cloud.nex
# 10 lines header
# + num taxa
# + line
# tress start
# don't forget to copy end


def parseNumber(string,word):
	string=str(string)
	findWord = string.find(word)
	number = string[:findWord]
	return float(number)

def subSetTreeFile(suffix,cloud1,cloud2,ratioName):
	for t in glob.glob(suffix):
		# Split off name, dependant on what your iterating handle files are
		listName = t.split("_")
		# Replace name with 500, remove existing file, open new empty file
		newName = listName[0]+"_"+str(ratioName)+"trees_"+listName[2]+"_"+listName[3]+"_"+listName[4]+"_"+listName[5]
		os.system("rm %s" % newName)
		outFile = open(newName,'a')
		# Get number of tips from file name
		numTips = parseNumber(listName[0],"tip")
		# Count lines to nexus header plus new size of first cloud
		firstCloud = 11 + int(numTips) + cloud1
		# Count lines to second cloud and size of new second cloud
		originalSize = parseNumber(listName[1],"trees")
		secondStart = 11 + int(numTips) + originalSize
		secondCloud = secondStart + cloud2
		# Iterate through file 
		print(newName)
		with open(t, 'r') as nexus:
			for i, line in enumerate(nexus):
				if i < firstCloud:
					#print(line)
					outFile.write(line)
				if secondCloud > i >= secondStart:
					#print(line)
					outFile.write(line)
			outFile.write("End;")
		outFile.close()

cloud1=1
cloud2=999
suffix="*cloud.nex"
ratioName="1v1000"
subSetTreeFile(suffix,cloud1,cloud2,ratioName)