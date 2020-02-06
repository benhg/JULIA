import os
import math
import time
import fnmatch
import datetime
from datetime import datetime 
import csv


def generate_data(proteomefile, directory):
    #code from Julia: cleans up files and creates individual files for each sequence

    with open(proteomefile) as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '>':
                line = line.replace(".", "")
                line = line.replace("\n","")
                line = line.replace("|", "_")
                line = line.split(" ")[0]
                line = line[1::]
                with open(f'{directory}/filenames.txt', 'a') as f: #this makes a new files called filenames.txt
                    f.write(line + '.fasta' + '\n')
            else :
                next
    with open(proteomefile) as f:
        lines = f.readlines()
        subfile = 0
        for line in lines:
            if line[0] == '>':
                line = line.replace(".", "")
                line = line.replace("|", "_") # must remove | and . characters, they confuse the computer
                line = line.split(" ")[0] # gets rid of annotations at the end of the seq name
                line = line[1::]
                line = line.replace("\n","")
                subfile = line + '.fasta'
                with open(f'{directory}/{subfile}', 'w') as f:
                    f.write('>' + line + '\n') 
            else :
                    with open(f'{directory}/{subfile}', 'a') as f:
                        f.write(line)
