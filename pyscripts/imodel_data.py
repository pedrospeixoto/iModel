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

import matplotlib

#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D


class imodelData(object):
	
	datafile = "errors.txt" #Data file name
	fancynames = {} #Naming convention
	options = {} #Graph structuring options
	varoptions = {} #For each variable option, the possible labels 
	filters = {} #For each variable from varoptions, a list of selected labels to use (or all)
	
	outloopname="OutLoop"
	inloopname="InLoop"
	midloopname="MidLoop"
	xvarname = "xVar"

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
				y = [ i.rstrip(string.digits) for i in y]
				#y = [ i.rstrip(string.digits).replace('_', '') for i in y]
			if xisnumber: 
				if all( i==int(i) for i in y): #Clean ints
					y = [ int(i) for i in y]
				else: #ignore floats
					continue
			
			self.data[x] = y

			#Sort and make unique to create list of options
			y=sorted(set(y))
			print(x, y)
			d[x]=y

		print()
		self.varoptions=d

		#Get all combinations for outer loops
		outloop = []
		outlabel = []
		for xout in self.options[self.outloopname]:
			outlabel.append(xout)
			outloop.append(self.varoptions[xout])

		self.outlabel=outlabel
		self.outopt=list(itertools.product(*outloop))
		print("Figures Options:")
		print(self.outlabel)
		print(self.outopt)
		print()

		#Get all combinations for inner  loops
		inloop = []
		inlabel = []
		for xout in self.options[self.inloopname]:
			inlabel.append(xout)
			inloop.append(self.varoptions[xout])
		self.inlabel=inlabel
		self.inopt=list(itertools.product(*inloop))
		print("In graph options:")
		print(self.inlabel)
		print(self.inopt)

		return

	def ConfigFigures(self):
		print("")
		print("Printing figures")
		print(self.options[self.midloopname][0])
		print(self.outlabel)
		
		for out in self.outopt:
			#Set title
			title=""
			for i, x in enumerate(out):
				y=self.fancynames.get(x,x)
				title=title+self.outlabel[i]+str(y)+" "
			print(title)

			#Filter data frame for this case
			datalocal=self.data
			for i, col in enumerate(out):
				print(self.outlabel[i],col)
				datalocal=datalocal.loc[datalocal[self.outlabel[i]] == col]
			print(datalocal)
			
			n=len(self.options[self.midloopname])
			figure = PlotterPanel( n, title, [self.options[self.xvarname]]*n, self.options[self.midloopname])
			c = 0
			#https://www.digitalocean.com/community/tutorials/data-analysis-and-visualization-with-pandas-and-jupyter-notebook-in-python-3
			
			for i, opt in enumerate(self.inopt): #Different Lines in graphs
				name=""
				for o in opt: #Join labels to get a full name
					print(o)
					
					name=name+self.fancynames.get(o, o).strip()+"_"
				name=name[0:len(name)-1]
				print(name)
				
				
				for j, pan in enumerate(self.options[self.midloopname]): #Panel 
					
					datapivot=pd.pivot_table(datalocal,  pan,  self.options[self.xvarname][0], self.inlabel)
					print(datapivot)	
					
					dataindex = datalocal.set_index(self.inlabel).sort_index()
					
					print(dataindex)
					print(dataindex.loc[self.inopt[0]][pan].values)

					sys.exit(1)		
					x = []
					y = []
					print()
					print(pan)
					for k, val in enumerate(self.data[pan]):
						include = 0
						for l, o in enumerate(opt):
							print(l,o,self.data[self.inlabel[l]][k]) 
							if o == self.data[self.inlabel[l]][k]:
								include = include + 1
						if include == len(opt):
							xtmp=self.data[self.options[self.xvarname][0]][k]
							print(xtmp,val)
							x.append(xtmp)
							y.append(val)
						
					figure.plot( j, x, y, label=name, i=c)
					c=c+1

			outname=self.datafile.replace('.txt','')
			name=name.replace(" ", "_")
			outname=outname+name+".eps"
							
			figure.finish(outname)
			plt.show()
	
class PlotterPanel(object):
	fontsize=16
	fontlarge=20
	fontsmall=14
	dimx=18
	dimy=8
	colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	linestyles = ['-', '--', ':', '-.', '-']
	markers = ['.', ',', 'o', 'v', '+', 'x']
	n=1
	
	def __init__(self, n, title, xlabel, ylabel):
		self.n=n
		self.fig, self.ax = plt.subplots(1, self.n, figsize=(self.dimx, self.dimy))
		for i in range(n):
			self.ax[i].set_xscale("log", nonposx='clip')
			self.ax[i].set_yscale("log", nonposy='clip')
			self.ax[i].set_title(title, fontsize=self.fontsize)
			self.ax[i].set_xlabel(xlabel[i], fontsize=self.fontsize)
			self.ax[i].set_ylabel(ylabel[i], fontsize=self.fontsize)
		
		#self.ax.set_xticks()
		#self.ax.set_yticks()

	def plot(self, pan, x, y, label, i):
		
		if len(x) == 0:
			return
		x, y = (list(t) for t in zip(*sorted(zip(x,y))))
		i=i % 5
		
		if pan > 0:
			print(label)
			self.ax[pan].plot(x, y, marker=self.markers[i], linestyle=self.linestyles[i], label=label)
		else:
			self.ax[pan].plot(x, y, marker=self.markers[i], linestyle=self.linestyles[i])
		
		return
	
	def finish(self, outname):
		self.fig.legend(fontsize=self.fontsmall, loc = "center right", borderaxespad=0.0)
		self.fig.subplots_adjust(right=0.85)
		self.fig.subplots_adjust(left=0.05)
		#for i in range(self.n):
		#	self.ax[i].set_frame_on(False)
		print(outname)
		self.fig.savefig(outname)
		return
	