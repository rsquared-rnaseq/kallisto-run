from os.path import join

configfile: "config.yaml"

# hardcode trials for now
# TRIAL="12-C-0114"
TRIAL="11-C-0123"
# TRIAL="14-C-0022"
SAMPLES, TUMORS = glob_wildcards(config["data_dir"].format(trial=TRIAL) + "/{samp}/Raw_fastq/{tumors}_1.fastq.gz")

# Directories
RAW_RESPONSE_DIR = join(config["data_dir"], "response_data")

# Files
RAW_METADATA_FILE = join(RAW_RESPONSE_DIR, "{trial}-tumor-response-sex-age.tsv")
PROCESSED_METADATA_FILE = join(config["data_dir"], "metadata.tsv")


rule all:
    input:
        expand(config["output_dir"].format(trial=TRIAL) + "/{samp}/{tumors}/abundance.tsv", zip, samp=SAMPLES, tumors=TUMORS)

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

rule generate_metadata_file:
    input:
        RAW_METADATA_FILE,
        expand(config["output_dir"].format(trial=TRIAL) + "/{samp}/{tumors}/abundance.tsv", zip, samp=SAMPLES, tumors=TUMORS, trial="{trial}")
    output:
        METADATA_FILE
    script:
        "src/create-metadata.py"
