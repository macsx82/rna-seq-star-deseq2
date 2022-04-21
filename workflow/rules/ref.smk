rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    resources:
        mem_mb=2000,
    wrapper:
        "v1.3.2/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/get_annotation.log",
    resources:
        mem_mb=2000,
    wrapper:
        "v1.3.2/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    resources:
        mem_mb=20000,
    envmodules:
        "samtools"
    wrapper:
        "v1.3.2/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        idx=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    params:
        algorithm="bwtsw",
    resources:
        mem_mb=40000,
    cache: True
    envmodules:
        "bwa"
    wrapper:
        "v1.3.2/bio/bwa/index"


rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 4
    params:
        extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    resources:
        mem_mb=20000,
    envmodules:
        "star"
    wrapper:
        "v1.3.2/bio/star/index"
