rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        ref_fasta=rules.get_genome.output[0],
        annotation=rules.get_annotation.output[0],
    output:
        config.get("paths").get("base_out") + "/1.STAR_ALIGN/pe/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        config.get("paths").get("base_out") + "/1.STAR_ALIGN/pe/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        config.get("paths").get("base_out") + "/1.STAR_ALIGN/pe/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-pe/{sample}-{unit}.log",
        "logs/star-pe/{sample}-{unit}.err"
    params:
        # index=config["star"]["star-genome"],
        extra=config["star"]["extra"],
        star_options=config["star"]["params"],
        tmpdir=config["paths"]["tmp"],
        out_prefix=config.get("star").get("base_out") + "/1.STAR_ALIGN/pe/{sample}-{unit}/"
    threads: 12
    resources:
        mem_mb=20000,
    envmodules:
        "star"
    shell:
        """
        temp=$(mktemp -u -d -p {params.tmpdir})
        STAR --runThreadN {threads} --genomeDir {input.ref_fasta} --sjdbGTFfile {input.annotation} --readFilesIn {input.fq1} {input.fq2} {params.star_options} {params.extra} --outTmpDir ${{temp}} --outFileNamePrefix {params.out_prefix} --outStd Log 1> {log[0]} 2> {log[1]}
        """
    # wrapper:
    #     "v1.3.2/bio/star/align"
        # "STAR "
        # " --runThreadN {snakemake.threads}"
        # " --genomeDir {index}"
        # " --readFilesIn {input_str}"
        # " {readcmd}"
        # " {extra}"
        # " "
        # " --outFileNamePrefix {tmpdir}/"
        # " --outStd Log "
        # " {log}"



# rule align_se:
#     input:
#         fq1=get_map_reads_input_R1,
#     output:
#         "results/star/se/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
#         "results/star/se/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
#         "results/star/se/{sample}-{unit}/ReadsPerGene.out.tab",
#     log:
#         "logs/star-se/{sample}-{unit}.log",
#     params:
#         index=config["star"]["star-genome"],
#         extra="--quantMode GeneCounts TranscriptomeSAM "
#         "--outSAMtype BAM SortedByCoordinate "
#         "--outFilterIntronMotifs RemoveNoncanonical "
#         "--chimSegmentMin 10 "
#         "--chimOutType SeparateSAMold "
#         "--outSAMunmapped Within "
#         "--sjdbGTFfile {} {}".format(
#             config["star"]["gtf"], config["star"]["params"]
#         ),
#     threads: 16
#     resources:
#         mem_mb=20000,
#     envmodules:
#         "star"
#     wrapper:
#         "v0.75.0/bio/star/align"

rule index_coord:
    input:
        get_star_bam,
    output:
        config.get("paths").get("base_out") + "/1.STAR_ALIGN/{ends}/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/{sample}-{unit}.{ends}.sortedByCoord.log"
    threads: 1        
    resources:
        mem_mb=20000,
    envmodules:
        "samtools"
    shell:
        """
        samtools index {input[0]} 2>&1 {log}
        """
