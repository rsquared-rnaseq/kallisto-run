configfile: "config.yaml"

import glob
import os

# Get all the FASTQs we need to process
fastqs_pair1 = glob.glob(config["data_dir"] + "/**/*_1.fastq.gz", recursive=True)
fastqs_pair2 = glob.glob(config["data_dir"] + "/**/*_2.fastq.gz", recursive=True)
if len(fastqs_pair1) != len(fastqs_pair2):
    print("Error: This pipeline only supports paired-end reads")
    quit()

fastqs = list(zip(sorted(fastqs_pair1), sorted(fastqs_pair2)))
outputs = []

for fq1, fq2 in fastqs:
    basename = os.path.basename(fq1)
    # basename[:4] returns the 4-digit tumor ID that the fastq belongs to.
    # basename[:-10] returns the filename without the extension
    outputs.append(config["output_dir"] + "/{tumor}/{outdir}".format(
        tumor=basename[:4],
        outdir=basename[:-len("_1.fastq.gz")]))
output_fns = [odir + "/abundance.tsv" for odir in outputs]

rule all:
    input: output_fns

rule run_kallisto:
    input: fqs=fastqs, index=config['index_file']
    output: output_fns
    run:
        for (fq1, fq2), odir in zip(fastqs, outputs):
            kallisto_cmd = "kallisto quant -i {index} -o {od} -b 100 {f1} {f2}".format(
                index=input.index,
                od=odir,
                f1=fq1,
                f2=fq2
            )
            shell(kallisto_cmd)