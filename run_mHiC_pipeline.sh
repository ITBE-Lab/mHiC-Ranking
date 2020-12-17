# create directory for slurm logs, pipeline will hang if it doesn't exist
# (the folder name is defined in cluster.yaml)
mkdir -p .slurm_logs

snakemake -j 100 --use-conda --cluster-config cluster.yaml --cluster \
"sbatch --partition={cluster.queue} \
--cpus-per-task={cluster.threads} \
--mem={cluster.memory} \
--output={cluster.output} \
--job-name={cluster.name} \
--time={cluster.time}"
