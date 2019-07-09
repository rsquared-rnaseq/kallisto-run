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
        expand(config["output_dir"] + "/{samp}/{tumors}/abundance.tsv", zip, samp=SAMPLES, tumors=TUMORS, trial=TRIALS)
rule dl_ref_transcriptome:
    output:
        ref_tscriptome=config["base_dir"] + "/ref_transcriptome.fq.gz"
    shell:
        "wget {config[ref_transcriptome_url]} -O {output.ref_tscriptome}"

rule create_kallisto_index:
    input:
        ref_tscriptome=config["base_dir"] + "/ref_transcriptome.fq.gz"
    output:
        index_file=config["base_dir"] + "/GRCh38.idx"
    shell:
        "kallisto index -i {output.index_file} {input.ref_tscriptome}"

rule move_fastq_files:
    input:
        fq1=config["data_dir"] + "/{samp}/{tumors}_1.fastq.gz",
        fq2=config["data_dir"] + "/{samp}/{tumors}_2.fastq.gz",
    output:
        fq1=config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_1.fastq.gz",
        fq2=config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_2.fastq.gz",
    shell:
        "mv {input.fq1} {output.fq1}; mv {input.fq2} {output.fq2}"

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
