import os
import parsl
from parsl.app.app import python_app, bash_app
from parsl.config import Config
from parsl.executors.threads import ThreadPoolExecutor

config = Config(executors=[ThreadPoolExecutor(max_threads=20)])

parsl.load(config)


@bash_app
def star_index(directory, genomeDir, filename):
  filename = filename.strip()
  indexingstar = f'STAR --runThreadN 1 --runMode genomeGenerate --genomeDir  "{genomeDir}" --genomeFastaFiles "{directory}{filename}" --genomeSAindexNbases 2'
  indexingstar = f"'{indexingstar}'"
  # indexing = f"SGE_Batch -r '{genomeDir}/run_dir' -c {indexingstar} -P1"
  print(indexingstar)
  return indexingstar
