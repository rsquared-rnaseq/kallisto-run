#!/bin/bash
module load snakemake kallisto || exit 1

snakemake --cluster "sbatch --time=08:00:00 --mem=16g --cpus-per-task=4" --jobs 100 --latency-wait 60 --keep-going --local-cores 4 all
