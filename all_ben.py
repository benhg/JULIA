import os
import sys
import math
import time
import fnmatch

import parsl
from parsl.app.app import python_app, bash_app
from parsl.config import Config
from parsl.providers import GridEngineProvider
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_route

from data_generation import generate_data

config = Config(
    executors=[HighThroughputExecutor(worker_debug=True,cores_per_worker=16,
                                      address=address_by_route(), 
                                      provider=GridEngineProvider( walltime='10000:00:00',nodes_per_block=1,
                                                                  init_blocks=1, 
                                                                   max_blocks=4, scheduler_options="#$ -pe smp 48"), 
                                      label="workers")
              ],
)


parsl.set_stream_logger() 
parsl.load(config)

proteomefile = sys.argv[1]
directory = f'/home/users/ellenrichards/{sys.argv[2]}/'
threshold = 1000

if not os.path.isdir(directory):
    os.makedirs(directory)

generate_data(proteomefile, directory)

os.system(f"touch {directory}/hopper.txt")
os.system(f"touch {directory}/final.txt")
os.system(f"touch {directory}/never.txt")
os.system(f"touch {directory}/maybehopper.txt")
os.system(f"touch {directory}/index_hopping_output.txt")
os.system(f"echo 'filename  uniquely  multi  totalReads  uniquelyAGAINST  multiAGAINST  totalreadsAGAINST  percRatio'  >> {directory}/index_hopping_output.txt")

@bash_app
def run_single_index(filename, directory):
    import os     
    threshold=1000
    filename = filename.strip()
    os.chdir(directory)
    mvalue = str(f"{filename}")
    mvalue = mvalue.split("/n")[0]
    mvalue= mvalue.split('.fasta')[0]

    svalue = filename.split("_")[0]
    svalue = str(svalue.split('s')[1])
    svalue  = int(svalue)
    strsvalue = str(svalue)

    if not os.path.isdir(f"{directory}/{mvalue}"):
        os.makedirs(f"{directory}/{mvalue}")
    genomeDir = f"{directory}{mvalue}/gd"
    if not os.path.isdir(genomeDir):
        os.makedirs(genomeDir)

    os.chdir(directory + mvalue)  #change directory to mvalue folder we just made

    indexingstar = f'STAR --runThreadN 1 --runMode genomeGenerate --genomeDir  "{genomeDir}" --genomeFastaFiles "{directory}{filename}" --genomeSAindexNbases 2'
    return indexingstar

                
def run_single_align(filename):    
    star_align(directory, mvalue, strsvalue, genomeDir)
    
    output = 0

    uniquely = os.popen("awk 'FNR == 9{print $6}' " + strsvalue + "Log.final.out").read()
    if uniquely == 0 : 
        uniquely = '1'
    uniquely = int(uniquely)

    multi = os.popen("awk 'FNR == 24{print $9}' " + strsvalue + "Log.final.out").read()
    if multi == '0': 
        multi = '1'
    multi = int(multi)

    totalreads = os.popen("awk 'FNR == 6{print $6}' " + strsvalue+ "Log.final.out").read()

    if uniquely + multi >= threshold:
        os.system("echo " + filename + f">> {directory}/final.txt") 
        output = 1
        uniquelyC = 'NO'
        multiC = 'NO'
        totalreadsC = 'NO'
        percRatio = 'NO'

    if svalue <12: 
        counter = 1
    else: 
        counter = 15

    while output  == 0 and counter > 0:
        os.system("echo" + filename + f" >> {directory}/maybehopper.txt")
        if counter == svalue: 
            counter +=1
        if counter == 20 or counter == 7: 
            counter = 0
        elif counter == 17: 
            counter = 19

        strcounter = str(counter)

        star_align(directory, mvalue, strcounter, genomeDir)

        threshold = 1000

        uniquelyC = os.popen("awk 'FNR == 9{print $6}' " + strcounter + "Log.final.out").read()
        if uniquelyC == '0': 
            uniquelyC = '1'
        uniquelyC = int(uniquelyC)

        multiC =  os.popen("awk 'FNR == 24{print $9}' " + strcounter + "Log.final.out").read()
        if multiC == '0': 
            multiC = '1'
        multiC = int(multiC)

        totalreadsC = os.popen("awk 'FNR == 6{print $6}' " + strcounter + "Log.final.out").read()

        percRatio = (multiC + uniquelyC)/(multi + uniquely)
        
        if (multiC + uniquelyC) > threshold and percRatio > 100:
            output = 1
            os.system("echo " + filename + f" >> {directory}/hopper.txt")
        else:
            if counter == 17: 
                counter = 19
            elif counter == 20 or counter == 7: 
                counter = 0
            else: 
                counter = counter + 1


    if output == 1: #it needs to be added to a file
        os.system("echo " + str(filename) + "  " +  str(uniquely) + "  " + str(multi) + "  " + str(totalreads) + "  " +  str(uniquelyC) + "  " +   str(multiC) + "  " + str(totalreadsC) + "  " +  str(percRatio)   + f" >> {directory}/index_hopping_output.txt")


    if output == 0:
        os.system("echo " + filename + f" >> {directory}/never.txt")

@bash_app
def star_align(filename, directory, inputs=[]):
    import os
    import fnmatch
    filename = filename.strip()
    os.chdir(directory)

    mvalue = str(f"{filename}").split('.fasta')[0]
    svalue = str(int(str(filename.split("_")[0].split('s')[1])))

    genomeDir = f"{directory}{mvalue}/gd"
    
    outfilenameprefix = directory + mvalue +"/"+ svalue

    if int(svalue) < 10:
        fullS = "s00" + str(svalue)
    else:
        fullS = "s0" + str(svalue)
    
    for rrfile in os.listdir('/home/users/ellenrichards/binfordlab/raw_reads/'):
        if fnmatch.fnmatch(rrfile, "*" + fullS + "*R1*.fastq"):
            rawread1 = rrfile
        if fnmatch.fnmatch(rrfile, "*" + fullS + "*R2*.fastq"):
            rawread2 = rrfile

    alignstar= f'STAR --runMode alignReads --runThreadN 16 --genomeDir "{genomeDir}"  --readFilesIn /home/users/ellenrichards/binfordlab/raw_reads/"{rawread1}"  /home/users/ellenrichards/binfordlab/raw_reads/"{rawread2}" --outFileNamePrefix "{outfilenameprefix}" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outTmpDir "{directory}/{mvalue}/align_tmp/"'

    return alignstar



if __name__ == "__main__":
    files = open(f"{directory}/filenames.txt")
    print("starting indexing")
    index_futures = []
    for file in files:
        index_futures.append(run_single_index(file, directory))
    print("done indexing")
    index_futures = [i.result() for i in index_futures]
    align_futures = []
    files = open(f"{directory}/filenames.txt")
    for index, file in enumerate(files):
        align_futures.append(star_align(file, directory, inputs=[index_futures[index]]))
    align_futures = [a.result() for a in align_futures]

