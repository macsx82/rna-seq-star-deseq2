#download the selected version of the reference genome using ensembl
rule get_genome:
    output:
        base_ref + "/" + ref_genome
        # "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    threads: 1
    resources:
        mem_mb=2000,
    wrapper:
        "v1.17.2/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        base_ref + "/" + ref_annots
        # "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    threads: 1
    log:
        "logs/get_annotation.log",
    resources:
        mem_mb=2000,
    wrapper:
        "v1.17.2/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        rules.get_genome.output[0]
        # "resources/genome.fasta",
    output:
        base_ref + "/" + ref_genome +".fai",
        # "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    threads: 1
    resources:
        mem_mb=20000,
    envmodules:
        "samtools"
    shell:
        """
        samtools faidx {input[0]} 2> {log[0]}
        """

rule bwa_index:
    input:
        rules.get_genome.output[0]
        # "resources/genome.fasta",
    output:
        idx=multiext(base_ref + "/" + ref_genome, ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
        # idx=multiext(base_ref + "/" + ref_genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    # params:
    #     algorithm="bwtsw",
    threads: 1
    resources:
        mem_mb=40000,
    cache: True
    envmodules:
        "bwamem2"
    wrapper:
        "v1.17.2/bio/bwa-mem2/index"


rule star_index:
    input:
        fasta=rules.get_genome.output[0],
        annotation=rules.get_annotation.output[0]
        # fasta="resources/genome.fasta",
        # annotation="resources/genome.gtf",
    output:
        directory(base_ref + "/star_genome"),
        touch(base_ref + "/star_genome/done.txt")
    threads: 8
    params:
        # extra=lambda input: " ".join("--sjdbGTFfile", input[1],"--sjdbOverhang 100"),
        extra="--sjdbGTFfile " + base_ref + "/" + ref_annots+ " --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    resources:
        mem_mb=20000,
    envmodules:
        "star"
    wrapper:
        "v1.17.2/bio/star/index"
