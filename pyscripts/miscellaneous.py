#-----------------------------------------------------------------------
# Some aux plotting routines and constants
# Luan Santos - october 2023
#-----------------------------------------------------------------------

import subprocess
import os.path
import matplotlib.colors as mcolors

#----------------------------------------------------------------------------
# data directory
datadir = "../data/"

# graphs directory
graphdir = "../graphs/"

# par directory
pardir = '../par/'

# imodel latlon grid size
nlat = 720
nlon = 1440

#----------------------------------------------------------------------------
# colorbar for clouds and precipitation
cmap_data = [(1.0, 1.0, 1.0),
             (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
             (0.0, 1.0, 1.0),
             (0.0, 0.8784313797950745, 0.501960813999176),
             (0.0, 0.7529411911964417, 0.0),
             (0.501960813999176, 0.8784313797950745, 0.0),
             (1.0, 1.0, 0.0),
             (1.0, 0.6274510025978088, 0.0),
             (1.0, 0.0, 0.0),
             (1.0, 0.125490203499794, 0.501960813999176),
             (0.9411764740943909, 0.250980406999588, 1.0),
             (0.501960813999176, 0.125490203499794, 1.0),
             (0.250980406999588, 0.250980406999588, 1.0),
             (0.125490203499794, 0.125490203499794, 0.501960813999176)]
colormap = mcolors.ListedColormap(cmap_data, 'precipitation')
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
def replace_line(filename, content, line_number):
    import re
    if os.path.exists(filename): # The file exists
        # Open the grid file
        file  = open(filename, "r")
        lines = file.readlines()

        # new content
        lines[line_number-1] = content+'\n'

        # Close the file
        file.close()

        # Write new file
        with open(filename, 'w') as file:
            for line in lines:
                file.write(line)

    else:   # The file does not exist
        print("ERROR in edit_file_line: file"+filename+" not found in /par.")
        exit()
#----------------------------------------------------------------------------

