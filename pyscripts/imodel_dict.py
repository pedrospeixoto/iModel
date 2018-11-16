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

class Dictionary(object):
	names = {}
	
	def __init__(self, filename):	
		with open(filename, mode='r') as infile:
			reader = csv.reader(infile)
			self.names = {rows[0]:rows[1] for rows in reader}
		print(self.names)