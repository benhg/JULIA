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


# Parsl config for use on LC's campus cluster
# You should be able to change this configuration and reproduce the same results
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

# Enable parsl logging if you want, but it prints out a lot of (useful) info
# parsl.set_stream_logger()
parsl.load(config)


def parsl_first_align(directory, spider_type):
    """
    Index and run first align, Submitting all tasks to Parsl executor at the beginning.
    
    This is a blocking call. It first performs all of the indexing (total of about 1 minute for 600 tasks)
    Next, it calls the alignment on each of those indexed transcripts (~12 min per alignment).
    
    :param directory - directory that stores output and input
    """
    files = [f.strip() for f in open(f"{directory}/filenames.txt").readlines()]

    # Start the alignment processes
    print("starting first alignment")
    align_futures = []
    for index, file in enumerate(files):
        align_futures.append(
            star_align(
                file, spider_type))
    # Wait for the alignment to finish
    align_futures = [a.result() for a in align_futures]
    print("First alignment finished")


@bash_app
def star_align(filename, spider_type):
    mvalue = str(f"{filename}").split('.fasta')[0]      
    svalue = str(int(str(filename.split("_")[0].split('s')[1])))
    alignstar = f'STAR --runMode alignReads --runThreadN 16 --genomeDir "{directory}{mvalue}/gd"  --readFilesIn /home/users/ellenrichards/binfordlab/raw_reads/"{spider_type}_R1.fastq"  /home/users/ellenrichards/binfordlab/raw_reads/"{spider_type}_R2" --outFileNamePrefix {directory}{mvalue}/{svalue} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outTmpDir "{directory}/{mvalue}/align_tmp/"'
    return alignstar



def setup():
    """
    Set up the run - 
        Get input arguments
        Create necessary directories
        Generate individual transcript FASTA files
        Create output files
    """
    directory = f'{sys.argv[1]}/'
    spider_type = f'{sys.argv[2]}'

    # Deletes output directory if it exists
    if not os.path.isdir(directory):
        print("Output directory must exist already, with filenames file")
        exit(1)

    main_output = f"{directory}/index_hopping_output.txt"

    Path(f"{directory}/hopper.txt").touch()
    Path(f"{directory}/final.txt").touch()
    Path(f"{directory}/never.txt").touch()
    Path(f"{directory}/maybehopper.txt").touch()
    Path(main_output).touch()

    with open(main_output, "a") as fh:
        fh.write('filename  uniquely  multi  totalReads  uniquelyAGAINST  multiAGAINST  totalreadsAGAINST  percRatio\n')

    return directory, spider_type


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: ./redundancy_alignment.py <output directory> <spider_type>")
        exit(1)

    directory, spider_type = setup()
    parsl_first_align(directory, spider_type)  # Blocking call waits on all first aligns
