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

#For data structuring
import itertools

class Figure(object):
		def __init__(self, param, names, infile):
			label=""
			title=""
			for i, xout in enumerate(names):
				label=label+"_"+xout+str(int(param[i]))
				title=title+" "+xout+" "+str(int(param[i]))
			self.label=label
			self.param=param
			self.names=names
			self.filename=
		
		def panels():


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
		self.GetOptions()

	def FancyNames(self, filename):	
		self.fancynames = {}
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			self.fancynames = {rows[0]:rows[1] for rows in reader}

	def GetOptions(self):
		#Check option in string variables
		print()
		print("String options to filter")
		d = {}
		for x in self.datastr:
			print(x)
			y=self.data[x]
			if x=='Grid': #Remove numbers from grid names
				y= [ i.rstrip(string.digits).replace('_', '') for i in y]
			y=sorted(set(y))
			print(y)
			d[x]=y
		print()
		print("Numerics options to filter")
		for x in self.datanum:
			
			y=self.data[x]
			y=sorted(set(y))
			if all( i==int(i) for i in y): #ignore floats
				print(x)
				d[x]=y
				print(y)	
			
		self.varoptions=d
		return

	def UserOptions(self, filename):	
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			#for row in reader:
			#	self.filter[row[0]]=row[1]
			userdata = list(reader)
			n=len(userdata)
			for i in range(n):
				if userdata[i][0] in self.varoptions.keys():
					#This option exists!
					if userdata[i][1] != "all":
						self.varoptions[userdata[i][0]]=userdata[i][1:]
				else:
					self.varoptions[userdata[i][0]]=userdata[i][1:]

	def ConfigFigures(self, input_filename):

		outloop = []
		for xout in self.varoptions["OutLoop"]:
			outloop.append(self.varoptions[xout])

		outeroptions=list(itertools.product(*outloop))
		#print(outeroptions)

		self.figures = []
		for op in outeroptions:
			self.figures.append(Figure(op, self.varoptions["OutLoop"]),input_filename)

		print(self.figures[1].label)
			
	

