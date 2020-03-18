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

from label_line import labelLinesEnd

class imodelData(object):
	
	fancynames = {} #Naming convention
	options = {} #Graph structuring options
	varoptions = {} #For each variable option, the possible labels 
	filters = {} #For each variable from varoptions, a list of selected labels to use (or all)
	
	yvarname="yVar"
	xvarname = "xVar"
	scalename = "Scale"

	def __init__(self, input_filenames): #Read data file
		self.datafiles = input_filenames
		self.data = []
		self.names = []
		print("Data file headers:")
		for file in self.datafiles:
			data = pd.read_csv(file, delim_whitespace=True)
			filename=os.path.splitext(file)[0]
			name=os.path.basename(filename)
			print(name)
			self.names.append(name)
			self.data.append(data)
			print(file, list(data))
		
		print()		

	def FancyNames(self, filename):	#Dictionary for fancy names
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			self.fancynames = {rows[0]:rows[1] for rows in reader}

		for i, name in enumerate(self.names):
			self.names[i]=self.fancynames.get(name,name)	

		print(self.names)

	def UserOptions(self, filename): #Plotting options
		# -yVar - defines the data to be plotted (yvar)
		# -xVar - variable for x-axis

		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			next(reader) #skip first line that has comments
			
			for row in reader:
				print(row[0], row[1:])
				self.options[row[0]]=row[1:]

			self.xvar=self.options[self.xvarname][0]
			self.logscale=[0,0]
			self.logscale[0]=self.options[self.scalename][0]=="log"
			self.logscale[1]=self.options[self.scalename][1]=="log"
			
			print()
			print("Graph options:")
			print(self.options)
			print()
	

	def BuildFigures(self, title):
		print("")
		print("Printing figures")
		print()
		
		
		n = len(self.options[self.yvarname])
		print("Number of panels: ", n)

		
		#Set pretty y-labels
		ylabels=self.options[self.yvarname].copy()
		for ilab, lab in enumerate(ylabels):
			ylabels[ilab]=self.fancynames.get(lab,lab)
		print(ylabels)
		
		figure = PlotterPanel( n, title, [self.xvar]*n, ylabels, self.logscale)

		#Plot data for each panel
		for i, pan in enumerate(self.options[self.yvarname]): #Panel 

			#Loop over data from different files and plot each line
			for j, data in enumerate(self.data):
				name=self.names[j]
				print(name, self.xvar, pan)
				
				x=data[self.xvar].values.T
				y=data[pan].values
				print(x,y)
				figure.plot(i, x, y, label=name, i=j)

			figure.ax[j].legend()

		title=title.replace(" ", "")
		outname=title+".eps"
						
		figure.finish(outname)

		plt.show()

		return figure

			
	
class PlotterPanel(object):
	fontsize=16
	fontlarge=20
	fontsmall=14
	dimx=20
	dimy=10
	colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	linestyles = ['-', '--', ':', '-.', '-']
	#markers = ['.', ',', 'o', 'v', '+', 'x']
	markers = ['', '', '', '', '', '']
	n=1
	
	def __init__(self, n, title, xlabel, ylabel, logscale):
		self.n=n
		#self.fig, self.ax = plt.subplots(1, self.n) #, figsize=(self.dimx, self.dimy))
		self.fig, self.ax = plt.subplots(1, self.n, figsize=(self.dimx, self.dimy))
		if n < 2:
			self.ax=[self.ax, self.ax] #This is just to allow a single panel figure with subscriptable axes

		for i in range(n):
			#log scale
			if logscale[0]:
				self.ax[i].set_xscale("log", nonposx='clip')
			else: #linear scale
				self.ax[i].set_xscale("linear", nonposx='clip')
			
			#log scale
			if logscale[1]:
				self.ax[i].set_yscale("log", nonposy='clip')
			else: #linear scale
				self.ax[i].set_yscale("linear", nonposy='clip')

			self.ax[i].set_title(title, fontsize=self.fontsize)
			self.ax[i].set_xlabel(xlabel[i], fontsize=self.fontsize)
			self.ax[i].set_ylabel(ylabel[i], fontsize=self.fontsize)
		

		#self.ax.set_xticks()
		#self.ax.set_yticks()

	def plot(self, pan, x, y, label, i):
		
		if len(x) == 0:
			return
		
		n=len(x)

		if pan > 0 or self.n==1:
			print(label)
			#n=6
			#self.ax[pan].plot(x, y, marker=self.markers[i], linestyle=self.linestyles[i], label=label)
			self.ax[pan].plot(x[0:n], y[0:n], marker=self.markers[i], linestyle=self.linestyles[i], label=label)
		else:
			#self.ax[pan].plot(x, y, marker=self.markers[i], linestyle=self.linestyles[i])
			self.ax[pan].plot(x[0:n], y[0:n], marker=self.markers[i], linestyle=self.linestyles[i])
		
		#Add label
		#dx=(x[n-1]-x[n-2])/5.0
		#xlab=x[n-1]+dx
		#ylab=y[n-1]
		#self.ax[pan].set_xlim(right=xlab+5*dx)
		#self.ax[pan].text(xlab,ylab,label,backgroundcolor=self.ax[pan].get_facecolor())

		return
	
	def finish(self, outname):
		#self.fig.legend(loc='upper left', bbox_to_anchor=(0.2, 0.90))
		#self.fig.legend(loc='best')
		#self.fig.subplots_adjust(right=0.85)
		#self.fig.subplots_adjust(left=0.05)
		#for i in range(self.n):
		#	self.ax[i].set_frame_on(False)
		self.fig.tight_layout()

		print(outname)
		self.fig.savefig(outname)
		return
	
