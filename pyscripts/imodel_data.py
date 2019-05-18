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

import pandas as pd

class imodelData(object):
	
	datafile = "errors.txt" #Data file name
	fancynames = {} #Naming convention
	options = {} #Graph structuring options
	varoptions = {} #For each variable option, the possible labels 
	filters = {} #For each variable from varoptions, a list of selected labels to use (or all)
	
	outloopname="OutLoop"
	inloopname="InLoop"
	midloopname="MidLoop"

	def __init__(self, input_filename): #Read data file
		self.datafile = input_filename
		self.data = pd.read_csv(input_filename, delim_whitespace=True)
		#print(self.data.dtypes)
		#print(self.data.values)
		print("Data file header:")
		print(list(self.data))
		print()		

	def FancyNames(self, filename):	#Dictionary for fancy names
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			self.fancynames = {rows[0]:rows[1] for rows in reader}

	def UserFilters(self, filename): #User defined filters
		#For each variable in header a filter can be defined selecting specific values
		self.filters = {}
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			next(reader) #skip first line that has comments
			
			for row in reader:
				self.filters[row[0]]=row[1:]
			print("Filters:")
			print(self.filters)
			print()

	def UserOptions(self, filename): #PLotting options
		# -Inloop - goes into a graph
		# -MidLoop - goes into separate panels - defines the flot data to be plotted
		# -OutLoop - goes into different figures
		# -xVar - variable for x-axis

		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			next(reader) #skip first line that has comments
			
			for row in reader:
				self.options[row[0]]=row[1:]
			print("Graph nesting:")
			print(self.options)
			print()
	
	def OrganizeOptions(self):
		#Get unique values from data to be used as options
		print()
		print("Options to filter")
		
		d = {}
		
		for x in self.data:

			y = self.data[x].values

			#Clean labels
			xisnumber=np.issubdtype(self.data[x].dtype, np.number)
			if not xisnumber: #Clean string names (strip digits)
				y = [ i.rstrip(string.digits).replace('_', '') for i in y]
			if xisnumber: 
				if all( i==int(i) for i in y): #Clean ints
					y = [ int(i) for i in y]
				else: #ignore floats
					continue
			
			#Sort and make unique
			y=sorted(set(y))
			print(x, y)
			d[x]=y

		print()
		self.varoptions=d

		#Get all combinations for inner and out loops
		outloop = []
		outlabel = []
		for xout in self.options[self.outloopname]:
			outlabel.append(xout)
			outloop.append(self.varoptions[xout])

		outopt=list(itertools.product(*outloop))
		print(outlabel)
		print(outopt)
		return

	def ConfigFigures(self):
		
		#Get all options for outer loop (multiple combinations of choices)
		outloop = []
		for xout in self.varoptions["OutLoop"]:
			outloop.append(self.varoptions[xout])

		outeroptions=list(itertools.product(*outloop))
		#print(outeroptions)

		self.figures = []
		for op in outeroptions:
			self.figures.append(Figure(op, self.varoptions["OutLoop"], self))
		print('Created figure layouts')
		for fig in self.figures:
			print(fig.param)
			for yvar in self.varoptions["MidLoop"]:
				print(yvar)
				fig.addpanel(yvar, self)
			
	
class Figure(object):
	panels = []
	def __init__(self, param, names, data):
		label=""
		title=""
		for i, xout in enumerate(names):
			label=label+xout+str(int(param[i]))+"_"
			title=title+xout+" "+str(int(param[i]))+" "
		self.label=label
		self.param = {}
		for i, name in enumerate(names):
			self.param[name]=param[i]
		self.filename=data.infile

	def addpanel(self, yvar, data):
		self.panels.append(Panel(yvar, self.param, data))

class Panel(object):
	def __init__(self, yvar, param, data):
		self.xvar=data.varoptions['xVar'][0]
		self.yvar=yvar
		print(self.xvar, self.yvar)
		self.x = []
		self.y = []
		for i, val in enumerate(data.data[yvar]):
			print(i)
			for figopt in param.keys():
				print(figopt)
				if data.data[figopt][i] != param[figopt]:
					print('skip this line', data.data[figopt][i])
					break
			print('added this line' )
			self.x.append(data.data[self.xvar][i])
			self.y.append(val)
		print(self.x,self.y)
			#or varop in data.varoptions['InLoop']:
			#	print(varop, fig.param.keys())
			#	if varop in fig.param.keys():
			#		print(varop, data.data[varop][i], fig.param[varop])
			#		if data.data[varop][i] != fig.param[varop] :
			#			print(i, val, data.data[self.xvar][i], val)
				#if field[i] == f and methods[i] == mtd and gridnames[i] == grd:
				#	x.append(gridres[i])
				#		ymax.append(maxerrors[i])
				#
				#figure.plot( 0, x, ymax, label=name, i=c)
				#figure.plot( 1, x, yrms, label=name, i=c)
				#c = c + 1
				
				#plt.show()