#!/bin/bash
#SBATCH --partition=slim18	# (queue, see sinfo)
#SBATCH --mem=4G 				#memory pool for all cores
#SBATCH --cpus-per-task=1                     # number of threads
#SBATCH --time=4-00:00:00

module load Anaconda3/2019.03
source activate snakemake
module load ngs/samtools/1.9 \
            ngs/bwa/0.7.16 \
            python/3.6

source /work/project/ladsie_002/myPython-3/bin/activate

cd /work/project/ladsie_002/analyses/scripts/2019-09-05_mHiC_pipeline

snakemake -j 30 --rerun-incomplete --cluster-config cluster.yaml --cluster \
"sbatch --partition={cluster.queue} \
--cpus-per-task={cluster.threads} \
--mem={cluster.memory} \
--output={cluster.output} \
--job-name={cluster.name} \
--time={cluster.time}"
