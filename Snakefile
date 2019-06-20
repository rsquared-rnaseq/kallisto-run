configfile: "config.yaml"

SAMPLES, TUMORS = glob_wildcards(config["data_dir"] + "/{samp}/Raw_fastq/{tumors}_1.fastq.gz")

rule all:
    input:
        expand(config["output_dir"] + "/{samp}/{tumors}/abundance.tsv", zip, samp=SAMPLES, tumors=TUMORS)

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