import glob
import os
import pathlib
import concurrent.futures as cf
import functools as ft
import subprocess as sp


def kallisto_done_callback(result, fastq_output):
    if result.exception() is None:
        # Create a file that indicates processing is completed for that pair of FASTQs
        pathlib.Path(os.path.join(fastq_output, ".proc_completed")).touch()


# TODO: Parse these all with argparse
data_base = "/Users/restifo/NIH/data/paired"
fastq_identifier = "*.fq"
output_dir = "/Users/restifo/NIH/data/paired/output"
index_file = "/Users/restifo/NIH/data/GRCh38.idx"

firstend_suffix = "_1.fq"
secondend_suffix = "_2.fq"
max_jobs = 4

fastqs = []
for unpaired_fastq in glob.glob(os.path.join(data_base, fastq_identifier)):
    # Get the base name of each FASTQ file so we can split them up into pair 1 and pair 2
    # fastq_basename = os.path.basename(unpaired_fastq).split(firstend_suffix)
    fastq_basename = unpaired_fastq.split(firstend_suffix)
    if len(fastq_basename) > 1:
        fastqs.append(fastq_basename[0])

# Make sure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Create a thread pool to manage launching kallisto processes. This also allows us to attach a callback to the process
# ending, which we will use for keeping track of data processing.
pool = cf.ThreadPoolExecutor(max_workers=max_jobs)
# Loop through each fastq file
for fastq in fastqs:
    # Make the kallisto output directory
    fastq_out = os.path.join(output_dir, os.path.basename(fastq))
    os.makedirs(fastq_out, exist_ok=True)

    # Build kallisto command
    kallisto_cmd = "kallisto quant -i {0} -o {1} -b 100 {2} {3}".format(index_file, fastq_out,
                                                                        fastq + firstend_suffix,
                                                                        fastq + secondend_suffix)
    print(kallisto_cmd)

    # Actually run kallisto (aka kall it lol)
    kall = pool.submit(sp.call, kallisto_cmd, shell=True)

    # Once kallisto is done running, this callback will write a file indicating this. This allows us to resume
    # processing a failed run by simply checking if this file exists
    kall.add_done_callback(ft.partial(kallisto_done_callback, fastq_out))
