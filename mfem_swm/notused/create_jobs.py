#! /usr/bin/python

import sys
import time
import commands
import os
import math

# working directory
working_directory=os.path.abspath('.')

femswe_cpu_set=[1,2,4,8,16,24,32]
femswe_cell_set=['hex','cube']

for cell in femswe_cell_set:
	for num_cpus in femswe_cpu_set:

		job_description_string='femswe_'+cell+'_'+str(num_cpus)
		job_filepath=working_directory+'/'+job_description_string+'.cmd'
		output_filepath=working_directory+'/'+job_description_string+'.out'
		job_starting_script='compare_femswe_'+cell+'.sh'

		job_file_content="""#! /bin/bash

# output
#SBATCH -o """+output_filepath+"""
# working directory
#SBATCH -D """+working_directory+"""
# job description
#SBATCH -J """+job_description_string+"""
#SBATCH --get-user-env
#SBATCH --partition=wsm
#SBATCH --ntasks=1-1
#SBATCH --cpus-per-task="""+str(num_cpus)+"""
#SBATCH --mail-type=end
#SBATCH --mail-user=d.gutermuth@t-online.de
#SBATCH --export=NONE
#SBATCH --time=00:30:00

source /etc/profile.d/modules.sh

export OMP_NUM_THREADS="""+str(num_cpus*2)+"""

export KMP_AFFINITY="compact"

./"""+job_starting_script+""" 8
"""

		#print job_file_content

		print "Writing jobfile '"+job_filepath+"'"
		f=open(job_filepath, 'w')
		f.write(job_file_content)
		f.close()
