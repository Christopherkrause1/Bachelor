import numpy as np
import glob

#merging different txt.files into one:
##################################################
new_t_1 = 0   #constant important for the first loop in the code
t_final = []  #arrays which are later filled with data in the loop
t_extra = []

merge = open("merge_file.txt", "w+") #creates merge file in which the data is merged

for text in glob.glob("Felix_Daten/*.txt"):            #loop over all the data files (loop goes through files in alphabetic order)

    fin = open(text, "r")    #opens each data to cache it
    t_data = np.genfromtxt(text, usecols=(0), unpack=True) # get timestamps from textfile
    new_t_2 = t_data - t_data[0] # convert stamps to seconds, which then start from 0
    t_final = np.append(t_final, [new_t_2 + new_t_1]) #array with all the times with will be plotted


    ###################
    t_extra = np.append(t_extra, new_t_2[-1]+new_t_1) #times of the measured damage rate
    ###################

    new_t_1 += t_data[-1]-t_data[0] #closes time gap between the files
    data = fin.read()
    fin.close()
    fout = open("merge_file.txt", "a") # opens the merge_file to save the data in it
    fout.write("#----------- data file: " + text +" --------------\n") #creates notification at the start of each file
    fout.write(data) #writes the cached data in it
    fout.close()
##################################################
