#! /usr/bin/env python3
#
#
#   Given data (x,y), add line plot to graph
#
#
#   Pedro Peixoto <ppeixoto@usp.br>
#
#
#---------------------------------------------------

import os
import sys
import stat
import math
import csv
import numpy as np
import string
import re 
import json


class Names(object):
	names = {}
	def __init__(self, filename):	
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			self.names = {rows[0]:rows[1] for rows in reader}
		#print(self.names)
		
#This filtering needs refinement for detailed filtering, now it is all or just one
class Filter(object):
	filter = {}
	def __init__(self, filename):	
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			#for row in reader:
			#	self.filter[row[0]]=row[1]
			self.filter = {rows[0]:rows[1] for rows in reader}
		

	def select(self, list, variable):
		listtmp=[]
		choice=self.filter.get(variable, "all")
		if choice=="all":
			listtmp=list
		else:
			for item in list:
				if choice in item:
					listtmp.append(item)
		print(variable, choice, listtmp)
		return listtmp

class imodelData(object):
	def __init__(self, input_filename):
		#Get header
		lines = open(input_filename).readlines()
		self.datahead = lines[0].split()
		#print(datahead)
		self.ncol=len(self.datahead)
		
		#Check first line to get types
		line1 = lines[1].split()
		#print(line1)
		datatypes = []
		datastr = []
		datanum = []
		for i, x in enumerate(line1):
			#print(x)
			try:
				float(x)
				datatypes.append('f16')	
				datanum.append(self.datahead[i])
			except ValueError:
				datatypes.append('U60')
				datastr.append(self.datahead[i])
			
		#Get data
		self.data = np.genfromtxt(input_filename, skip_header=1, dtype=datatypes, autostrip=True, names=self.datahead)
		self.datastr=datastr
		self.datanum=datanum

		print("Data with strings")
		print(datastr)
		print()
		print("Data with numbers")
		print(datanum)
		print(	)

	def getoptions(self):
		#Check option in string variables
		print("Options to filter")
		d = {}
		for x in self.datastr:
			print(x)
			y=self.data[x]
			if x=='Grid': #Remove numbers from grid names
				y= [ i.rstrip(string.digits).replace('_', '') for i in y]
			y=sorted(set(y))
			print(y)
			d[x]=y
		self.varoptions=d
		return