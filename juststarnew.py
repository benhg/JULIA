import os
import math
import time
import fnmatch
import datetime
from datetime import datetime 
import csv

proteomefile = '80percnuc.fasta'
#proteomefile = 'try7.txt'
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
            with open('filenames.txt', 'a') as f: #this makes a new files called filenames.txt
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
            with open(subfile, 'w') as f:
                f.write('>' + line + '\n') 
        else :
                with open(subfile, 'a') as f:
                    f.write(line)


os.system("touch hopper.txt")
os.system("touch final.txt")
os.system("touch never.txt")
os.system("touch maybehopper.txt")
os.system("touch index_hopping_output.txt")

csv_filename = "index_hopping_output.csv"
csv_fh = open(csv_filename, "write")
csv_writer = csv.DictWriter(csv_fh, fieldnames=["filename", "uniquely", "multi", "totalReads", "uniquelyAGAINST", "multiAGAINST", "totalReadsAGAINST", "percratio"])
csv_row = {"filename": "filename ", "uniquely": "uniquely ", "multi": "multi ", "totalReads": "totalreads ", "uniquelyAGAINST": "uniquelyC ", "multiAGAINST": "multiC ", "totalReadsAGAINST": "totalreadsC ", "percratio": "percRatio "}
csv_writer.writerow(csv_row)
os.system("echo 'filename  uniquely  multi  totalReads  uniquelyAGAINST  multiAGAINST  totalreadsAGAINST  percRatio'  >> index_hopping_output.txt")

directory = '/home/users/ellenrichards/'

for filename in os.listdir(directory):
    os.chdir(directory)
    if filename.endswith(".fasta") and filename.startswith("s0"):
        #finds and sets "mvalue" from filename
        #mvalue = filename.split("m")[1]
        mvalue = filename
        mvalue = mvalue.split("/n")[0]
        #mvalue = mvalue.split(".")[1]
        mvalue = mvalue.split(".fasta")[0]
        mvalue = str(mvalue)

        #finds and sets svalue from filename0
        fullS = filename.split("_")[0]
        strsvalue = fullS.split("s")[1]
        svalue = int(strsvalue)
        strsvalue = str(svalue)

        #this creates the folder and subfolder for each file
        os.system("mkdir " + mvalue)
        os.system("mkdir " + mvalue + "/gd")

        #this changes the directory to INSIDE the folder
        os.path.join(directory + mvalue + "/gd")
        genomeDir = directory + mvalue + "/gd"
        os.chdir(directory + mvalue)

        #Indexing
        #indexingstar='\"STAR --runThreadN 1 --runMode genomeGenerate --genomeDir \"' + genomeDir + '\" --genomeFastaFiles \"' + foldername +'\"/\"' + filename + '\" --genomeSAindexNbases 2\"'
        indexingstar='\"STAR --runThreadN 1 --runMode genomeGenerate --genomeDir \"' + genomeDir + '\" --genomeFastaFiles \"' + directory + filename + '\" --genomeSAindexNbases 2\"'
        indexing = "SGE_Batch -r 'genome_dir' -c " + indexingstar + " -P1"
        os.system(indexing)
        
        #set raw read files: use files that start with the fullS value
        for rawreadfile in os.listdir('/home/users/ellenrichards/binfordlab/raw_reads/'):
            if fnmatch.fnmatch(rawreadfile, "*" + fullS + "*R1*.fastq"): 
                rawread1= rawreadfile
            if fnmatch.fnmatch(rawreadfile, "*" + fullS + "*R2*.fastq"): 
                rawread2= rawreadfile
        
        #wait to make sure indexing has completed before continuing
        while not os.path.exists("Log.out"): 
            time.sleep(5) 
            
        #STAR against self-- get ownperc (also output into file)
        outfilenameprefix = directory + mvalue +"/"+ strsvalue 
        alignownstar= '\"STAR --runMode alignReads --runThreadN 15 --genomeDir \"' + genomeDir + '\" --readFilesIn /home/users/ellenrichards/binfordlab/raw_reads/\"' + rawread1 + '\"  /home/users/ellenrichards/binfordlab/raw_reads/\"' + rawread2 + '\" --outFileNamePrefix \"' + outfilenameprefix +'\" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 \"'
        aligning = "SGE_Batch -r align" + strsvalue + " -c " + alignownstar + " -P 15"
        os.system(aligning)
                
        #wait to make sure alignment has completed before continuing
        while not os.path.exists(strsvalue + "Log.final.out" ): 
            time.sleep(15) 
           
        #collect output data: from __Log.final.out output file
        uniquely = os.popen("awk 'FNR == 9{print $6}' " + strsvalue + "Log.final.out").read()
        if uniquely == '0': uniquely = '1'
        else: uniquely = uniquely.strip()
        uniquely = int(uniquely)
        multi = os.popen("awk 'FNR == 24{print $9}' " + strsvalue + "Log.final.out").read()
        if multi == '0': multi = '1'
        else: multi = multi.strip()
        multi = int(multi)
        totalreads = os.popen("awk 'FNR == 6{print $6}' " + strsvalue+ "Log.final.out").read()
        if totalreads == '0': totalreads = '1'
        else: totalreads = totalreads.strip()
        totalreads = int(totalreads)
        totalhits = uniquely + multi
        if totalhits == 0: totalhits = 1

        print totalhits
        
        threshold = 1000 #threshold = min number of hits to be sure it is its own hit
        output = 0 #sets output at 0, has not found a value for outputting yet

        if (totalhits) >= threshold:  #if the output is high enough, this means it is NOT a hopper, belongs with this taxon source, adds to a file for that name. 
            os.system("echo " + filename + " >> ../final.txt")
            uniquelyC = "NO" 
            multiC = "NO" 
            totalreadsC = "NO" 
            percRatio = "NO"
            output = 1 #needs to be outputted
            
        #sets the first taxon to test against
        if svalue <12: counter = 1
        if (svalue >=12): counter = 15
        
        while (output != 1 and counter > 0):
            os.system("echo " + filename + ">> ../maybehopper.txt") #if entering this into this loop, could be a hopper.

            if counter == svalue: 
                counter = counter + 1 #if the counter is the same as the self that just ran, no need to run again
                if counter == 20 or counter == 7: counter = 0 #make sure that this new number is not accurate
                elif counter == 17: counter = 19

            strcounter=str(counter) #this will be helpful later, when inputting the counter into STAR

            #this is to set the fullC value (so that I can call rawread values)
            if counter>9: fullC= "s0" + strcounter
            elif counter<=9: fullC="s00" + strcounter
            
            #find raw read files
            for rawreadfile in os.listdir('/home/users/ellenrichards/binfordlab/raw_reads/'):
                if fnmatch.fnmatch(rawreadfile, "*" + fullC + "*R1*.fastq"): 
                    rawreadC1= rawreadfile
                if fnmatch.fnmatch(rawreadfile, "*" + fullC + "*R2*.fastq"): 
                    rawreadC2= rawreadfile

            #STAR against counter
            outfileCnameprefix = directory + mvalue +"/"+ strcounter 
            alignCstar= '\"STAR --runMode alignReads --runThreadN 15 --genomeDir \"' + genomeDir + '\" --readFilesIn /home/users/ellenrichards/binfordlab/raw_reads/\"' + rawreadC1 + '\"  /home/users/ellenrichards/binfordlab/raw_reads/\"' + rawreadC2 + '\" --outFileNamePrefix \"' + outfileCnameprefix +'\" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 \"'
            aligningC = "SGE_Batch -r align" + strcounter + " -c " + alignCstar + " -P 15"
            os.system(aligningC)
                
            #wait to make sure alignment has completed before continuing
            while not os.path.exists(strcounter + "Log.final.out" ): 
                time.sleep(15) 
                    
            #collect output data: from __Log.final.out output file
            uniquelyC = os.popen("awk 'FNR == 9{print $6}' " + strcounter + "Log.final.out").read()
            if uniquelyC == '0': uniquelyC = '1'
            else: uniquelyC = uniquelyC.strip()
            uniquelyC = int(uniquelyC)

            multiC =    os.popen("awk 'FNR == 24{print $9}' " + strcounter + "Log.final.out").read()
            if multiC == '0': multiC = '1'
            else: multiC = multiC.strip()
            multiC = int(multiC)

            totalreadsC = os.popen("awk 'FNR == 6{print $6}' " + strcounter + "Log.final.out").read()
            if totalreadsC == '0': totalreadsC = '1'
            else: totalreadsC = totalreadsC.strip()
            totalreadsC = int(totalreadsC)

            print uniquelyC
            print multiC

            totalhitsC = uniquelyC + multiC
            if totalhitsC == 0: totalhitsC = 1

            print totalhits
            print totalreads
            print totalhitsC
            print totalreadsC

            percRatio = ((totalhits/totalreads)*100.0)/((totalhitsC/totalreadsC)*100.0)
                                
            if (percRatio>=100 and (totalhitsC>=threshold)): #if it is a hopper, must be added to a HOPPER file
               output = 1 #stops loop
               os.system("echo " + str(filename) + "  " +  str(uniquely) + "  " + str(multi) + "  " + str(totalreads) + "  " +  str(uniquelyC) + "  " +   str(multiC) + "  " + str(totalreadsC) + "  " +  str(percRatio)   + "  >> ../hopper.txt")
                    
            else : #if it is NOT a hopper, continue and counter is added
                if counter == 17: counter = 19
                elif counter == 20 or counter == 7: counter = 0
                else : counter = counter + 1

        if output == 1: #it needs to be added to a file
            csv_row = {"filename": str(filename), "uniquely": str(uniquely), "multi": str(multi), "totalReads": str(totalreads), "uniquelyAGAINST": str(uniquelyC), "multiAGAINST": str(multiC), "totalReadsAGAINST": str(totalreadsC), "percratio": str(percRatio)}
            csv_writer.writerow(csv_row)
            os.system("echo " + str(filename) + "  " +  str(uniquely) + "  " + str(multi) + "  " + str(totalreads) + "  " +  str(uniquelyC) + "  " +   str(multiC) + "  " + str(totalreadsC) + "  " +  str(percRatio)   + "  >> ../index_hopping_output.txt")

        elif output==0: #was not found as a hopper and did not meet threshold for its own
            os.system("echo " + filename + " >> ../never.txt")
