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

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import re


class Plotter(object):
	fontsize=16
	dimx=10
	dimy=10
	colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	linestyles = ['-', '--', ':', '-.', '-']
	markers = ['.', ',', 'o', 'v', '+', 'x']
	
	def __init__(self, title, xlabel, ylabel):
		self.fig, self.ax = plt.subplots(figsize=(self.dimx, self.dimy))
		self.ax.set_xscale("log", nonposx='clip')
		self.ax.set_yscale("log", nonposy='clip')
		self.ax.set_title(title)
		self.ax.set_xlabel(xlabel, fontsize=self.fontsize)
		self.ax.set_ylabel(ylabel, fontsize=self.fontsize)
		#self.ax.set_xticks()
		#self.ax.set_yticks()

	def plot(self, x, y, label, i):
		print(label)
		if len(x) == 0:
			return
		x, y = (list(t) for t in zip(*sorted(zip(x,y))))
		i=i % 5
		
		self.ax.plot(x, y, marker=self.markers[i], linestyle=self.linestyles[i], label=label)
		
	def finish(self, outname):
		self.fig.legend(fontsize=15)
		print(outname)
		self.fig.savefig(outname)
		return
	
	
	
