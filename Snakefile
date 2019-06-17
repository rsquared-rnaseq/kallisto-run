DATA_DIR = "/data/Robinson-SB"
samples, subds, fqnames, endings = glob_wildcards(DATA_DIR + "/{sample}/{subd}/{fqname}_{ending}.fastq.gz")