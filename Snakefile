configfile: "config.yaml"

TRIALS1, SAMPLES1, TUMORS1 = glob_wildcards(config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_1.fastq.gz")
TRIALS2, SAMPLES2, TUMORS2 = glob_wildcards(config["data_dir"] + "/{samp}/{tumors}_1.fastq.gz")

TRIALS = TRIALS1+TRIALS2
SAMPLES= SAMPLES1+SAMPLES2
TUMORS = TUMORS1+TUMORS2

rule all:
    input:
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
        fq2=config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_2.fastq.gz",
        index_file=config["base_dir"] + "/GRCh38.idx"
    output:
        ofile=config["output_dir"] + "/{samp}/{tumors}/abundance.tsv"
    params:
        odir=config["output_dir"] + "/{samp}/{tumors}"
    shell:
        "kallisto quant -i {input.index_file} -o {params.odir} -t {config[threads_per_proc]} -b 100 {input.fq1} {input.fq2}"
