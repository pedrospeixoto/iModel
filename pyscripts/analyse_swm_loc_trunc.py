#! /usr/bin/env python3
#---------------------------------
#   Plots errors of experiments
#      obtained from iModel output
#   Pedro Peixoto (ppeixoto@usp.br)
#   Novembre 2018
#----------------------------------


import sys
import os
import re
import string
import numpy as np

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

#Custom plotting setup
import imodel_plot
from imodel_plot import Plotter, PlotterPanel

import imodel_dict
from imodel_dict import Names, Filter

# input filename
input_filename = 'errors.txt'
if len(sys.argv) <= 1 :
	print("I need 1 argument:")
	print("A filename containing the errors generated by imodel")
	sys.exit(1)
	
if len(sys.argv) > 1:
	input_filename = sys.argv[1]
	

with open(input_filename) as f:
	lines = f.readlines()

#get header
head = lines[0]
#print(head)
if head[-1] == '\n':
	head = head[0:-1]
head = head.split()
print(head)

imethod = head.index("Methods")
ioperator = head.index("Operator")
igridname = head.index("Grid")
igridres = head.index("Mesh")
imaxerror = head.index("MaxError")
imaxerrorrel = head.index("MaxErrorRel")
irmserror = head.index("RMSError")

methods = []
operators = []
maxerrors = []
rmserrors = []
gridres = []
gridnames = []

for i, l in enumerate(lines[1:]):
	if l[-1] == '\n': #get rid of \n
		l = l[0:-1]

	d = l.split() #line data
	#print(d)
	# skip invalid nan's
	if d[imaxerror] == 'nan':
		continue
		
	operators.append(d[ioperator])
	methods.append(d[imethod])
	maxerrors.append(float(d[imaxerror]))
	rmserrors.append(float(d[irmserror]))
	gridres.append(int(d[igridres]))
	gridnames.append((d[igridname].rstrip(string.digits)))


#Get unique operator names
operators_list=sorted(set(operators))
methods_list=sorted(set(methods))
grids_list=sorted(set(gridnames))
print("Full options available")
print(operators_list)
print(methods_list)
print(grids_list)
print("---------------------")
print()


dict=Names("naming_conv.csv") #Load naming convention
filter=Filter("filter.csv") #load filtering conditions
print(dict)
#Apply filters
operators_list=filter.select(operators_list, 'operator')
grids_list=filter.select(grids_list, 'grid')
methods_list=filter.select(methods_list, 'method')
#methods_list=filter.select(methods_list)
#grids_list=filter.select(grids_list)
print("Filtred lists")
print(operators_list)
print(methods_list)
print(grids_list)
print("---------------------")
print()


#Plot for each operator
for oper in operators_list:
	outname=input_filename.replace('.txt', "_"+oper+".eps")
	#outnamerms=input_filename.replace('.txt', "_"+oper+"_rms.eps")
	title=dict.names.get(oper, oper)
	figure = PlotterPanel( 2, title, ["grid points", "grid points"], ["max error", "rms error"])
	#figurerms = Plotter(oper, "grid points", "rms error")
	c = 0
	for mtd in methods_list:
		for grd in grids_list:
			name=mtd+"_"+grd[0:-1]
			name=dict.names.get(name, name)
			x = []
			ymax = []
			yrms = []
			for i, val in enumerate(maxerrors):
				if operators[i] == oper and methods[i] == mtd and gridnames[i] == grd:
					x.append(gridres[i])
					ymax.append(maxerrors[i])
					yrms.append(rmserrors[i])
			figure.plot( 0, x, ymax, label=name, i=c)
			figure.plot( 1, x, yrms, label=name, i=c)
			c = c + 1
			
			#plt.show()
			
	figure.finish(outname)

	#figurerms.finish(outnamerms)


plt.show()

