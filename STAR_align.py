import os
import fnmatch

def star_align(directory, mvalue, svalue, genomeDir):
    svalue = int(svalue)
    svalue = str(svalue)
    outfilenameprefix = directory + mvalue +"/"+ svalue

    if int(svalue) < 10: 
        fullS = "s00" + str(svalue)
    else:
        fullS = "s0" + str(svalue)
    print(fullS)
    for rrfile in os.listdir('/home/users/ellenrichards/binfordlab/raw_reads/'):
        if fnmatch.fnmatch(rrfile, "*" + fullS + "*R1*.fastq"): 
            rawread1= rrfile
        if fnmatch.fnmatch(rrfile, "*" + fullS + "*R2*.fastq"):
            rawread2= rrfile

    alignstar= f'STAR --runMode alignReads --runThreadN 15 --genomeDir "{genomeDir}"  --readFilesIn /home/users/ellenrichards/binfordlab/raw_reads/"{rawread1}"  /home/users/ellenrichards/binfordlab/raw_reads/"{rawread2}" --outFileNamePrefix "{outfilenameprefix}" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000'
    alignstar = f"'{alignstar}'"

    align = f"SGE_Batch -r {align_svalue} -c  {alignstar} -P 15"

    os.system(align)
