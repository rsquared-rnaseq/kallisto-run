configfile: "config.yaml"


# TRIAL="12-C-0114"
TRIAL="11-C-0123"
# TRIAL="14-C-0022"

SAMPLES, TUMORS = glob_wildcards(config["data_dir"].format(trial=TRIAL) + "/{samp}/Raw_fastq/{tumors}_1.fastq.gz")


rule all:
    input:
        expand(config["output_dir"].format(trial=TRIAL) + "/{samp}/{tumors}/abundance.tsv", zip, samp=SAMPLES, tumors=TUMORS)

# rule move_fastq_files:
#     input:
#         fq1=config["data_dir"] + "/{samp}/{tumors}_1.fastq.gz",
#         fq2=config["data_dir"] + "/{samp}/{tumors}_2.fastq.gz",
#     output:
#         fq1=config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_1.fastq.gz",
#         fq2=config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_2.fastq.gz",
#     shell:
#         "mv {input.fq1} {output.fq1}; mv {input.fq2} {output.fq2}"

rule process_fastq:
    input:
        fq1=config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_1.fastq.gz",
        fq2=config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_2.fastq.gz"
    output:
        ofile=config["output_dir"] + "/{samp}/{tumors}/abundance.tsv"
    params:
        odir=config["output_dir"] + "/{samp}/{tumors}"
    shell:
        "kallisto quant -i {config[index_file]} -o {params.odir} -t {config[threads_per_proc]} -b 100 {input.fq1} {input.fq2}"
