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

class Names(object):
	names = {}
	
	def __init__(self, filename):	
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			self.names = {rows[0]:rows[1] for rows in reader}
		print(self.names)
		
class Filter(object):
	filter = {}
	def __init__(self, filename):	
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			self.filter = {rows[0]:rows[1] for rows in reader}
	
	def select(self, list, option):
		listtmp=[]
		choice=self.filter.get(option, "all")
		if choice=="all":
			return list
		else:
			for item in list:
				if choice in item:
					listtmp.append(item)
			return listtmp
			