import numpy as np
import glob
import sys
import os.path

folder = input("Please enter the path in which the data is stored: ")   #input folder from the command line


#merging different txt.files into one: merge_file.txt
##################################################
new_t_1 = 0         #constant important for the first loop in the code
t_merged_data = []  #array which is later filled with data in the loop
t_extra = []        #array which is later filled with data in the loop

merge = open("merge_file.txt", "w+")                     #creates merge file in which the data is merged

if os.path.exists(folder) == True:                       #checks if the given folder exists
    for text in glob.glob(folder + "/*.txt"):            #loop over all the data files (loop goes through files in alphabetic order)

        fin = open(text, "r")                                         #opens each data to cache it
        t_data = np.genfromtxt(text, usecols=(0), unpack=True)        #get timestamps from textfile
        new_t_2 = t_data - t_data[0]                                  #convert stamps to seconds, which then start from 0
        t_merged_data = np.append(t_merged_data, [new_t_2 + new_t_1]) #array with all the times with will be plotted


        ###################
        #t_extra = np.append(t_extra, new_t_2[-1]+new_t_1) #times of the measured damage rate
        ###################

        new_t_1 += t_data[-1]-t_data[0]    #closes time gap between the files
        data = fin.read()
        fin.close()
        fout = open("merge_file.txt", "a") # opens the merge_file to save the data in it
        fout.write("#----------- data file: " + text +" --------------\n") #creates notification at the start of each file
        fout.write(data) #writes the cached data in it
        fout.close()
else:
    print('path not found') #Stops the program if the folder does not exist
    exit()
##################################################
