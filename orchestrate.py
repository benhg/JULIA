#!/usr/local/bin/env python

import os
import sys
from pathlib import Path
import shutil

import parsl
from parsl.app.app import bash_app
from parsl.config import Config
from parsl.providers import GridEngineProvider
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_route

from data_generation import generate_data

config = Config(
    executors=[HighThroughputExecutor(worker_debug=True,
                                      cores_per_worker=16,
                                      address=address_by_route(),
                                      provider=GridEngineProvider(walltime='10000:00:00',
                                                                  nodes_per_block=1,
                                                                  init_blocks=1,
                                                                  max_blocks=4,
                                                                  scheduler_options="#$ -pe smp 48"
                                                                  ),
                                      label="workers")
               ],
)


parsl.set_stream_logger()
parsl.load(config)


@bash_app
def run_single_index(filename, directory, threshold=1000):
    import os

    mvalue = str(f"{filename}").split("/n")[0].split('.fasta')[0]
    genomeDir = f"{directory}{mvalue}/gd"
    if not os.path.isdir(f"{directory}/{mvalue}"):
        os.makedirs(f"{directory}/{mvalue}")
    if not os.path.isdir(genomeDir):
        os.makedirs(genomeDir)

    # change directory to mvalue folder we just made
    os.chdir(f"{directory}/{mvalue}")

    indexingstar = f'STAR --runThreadN 1 --runMode genomeGenerate --genomeDir  "{genomeDir}" --genomeFastaFiles "{directory}{filename}" --genomeSAindexNbases 2'
    return indexingstar


@bash_app
def star_align(filename, directory, inputs=[]):
    import os
    import fnmatch
    filename = filename.strip()
    os.chdir(directory)

    mvalue = str(f"{filename}").split('.fasta')[0]
    svalue = str(int(str(filename.split("_")[0].split('s')[1])))

    genomeDir = f"{directory}{mvalue}/gd"

    outfilenameprefix = directory + mvalue + "/" + svalue

    if int(svalue) < 10:
        fullS = "s00" + str(svalue)
    else:
        fullS = "s0" + str(svalue)

    for rrfile in os.listdir(
            '/home/users/ellenrichards/binfordlab/raw_reads/'):
        if fnmatch.fnmatch(rrfile, "*" + fullS + "*R1*.fastq"):
            rawread1 = rrfile
        if fnmatch.fnmatch(rrfile, "*" + fullS + "*R2*.fastq"):
            rawread2 = rrfile

    alignstar = f'STAR --runMode alignReads --runThreadN 16 --genomeDir "{genomeDir}"  --readFilesIn /home/users/ellenrichards/binfordlab/raw_reads/"{rawread1}"  /home/users/ellenrichards/binfordlab/raw_reads/"{rawread2}" --outFileNamePrefix "{outfilenameprefix}" --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outTmpDir "{directory}/{mvalue}/align_tmp/"'

    return alignstar


def parsl_first_align(directory):
    files = [f.strip() for f in open(f"{directory}/filenames.txt").readlines()]

    # Start the indexing processes
    print("starting indexing")
    index_futures = []
    for file in files:
        index_futures.append(run_single_index(file, directory))
    # Wait for the indexing to finish
    index_futures = [i.result() for i in index_futures]
    print("done indexing")

    # Start the alignment processes
    align_futures = []
    for index, file in enumerate(files):
        align_futures.append(
            star_align(
                file, directory, inputs=[
                    index_futures[index]]))
    # Wait for the alignment to finish
    align_futures = [a.result() for a in align_futures]


def setup():
    proteomefile = sys.argv[1]
    directory = f'{sys.argv[2]}/'

    # Deletes output directory if it exists
    if os.path.isdir(directory):
        shutil.rmtree(directory, ignore_errors=True)
    
    os.makedirs(directory)

    generate_data(proteomefile, directory)

    main_output = f"{directory}/index_hopping_output.txt"

    Path(f"{directory}/hopper.txt").touch()
    Path(f"{directory}/final.txt").touch()
    Path(f"{directory}/never.txt").touch()
    Path(f"{directory}/maybehopper.txt").touch()
    Path(main_output).touch()

    with open(main_output, "a") as fh:
        fh.write('filename  uniquely  multi  totalReads  uniquelyAGAINST  multiAGAINST  totalreadsAGAINST  percRatio\n')

    return directory


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: ./orchestrate.py <path to input FASTA> <output directory>")
        exit(1)

    directory = setup()
    parsl_first_align(directory)  # Blocking call waits on all first aligns
