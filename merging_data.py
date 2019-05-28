import numpy as np
import glob, os

#merging different txt.files into one:
##################################################

new_t_1 = 0
merge = open("merge_file.txt", "w+") #creates merge file in which the data is merged
t_final = []
t_start = 1537644542.5356536

for text in glob.glob("Felix_Daten/*.txt"):

    fin = open(text, "r")    #opens each data to cache it
    # if 'new_t_1' in locals():
    # else:
    #new_t_1 = t_data[:-1]
    t_data = np.genfromtxt(text, usecols=(0), unpack=True) # get timestamps from textfile
    new_t_2 = t_data - t_data[0] # convert stamps to seconds, which then start from 0
    t_final = np.append(t_final, [new_t_2 + new_t_1])

    new_t_1 += t_data[-1]-t_data[0]
    data = fin.read()
    fin.close()
    fout = open("merge_file.txt", "a") # opens the merge_file to save the data in it
    fout.write("#----------- data file: " + text +" --------------\n") #creates notification at the start of each file
    fout.write(data) #writes the cached data in it
    fout.close()
##################################################
