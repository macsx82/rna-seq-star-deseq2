rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
    output:
        "results/star/pe/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/pe/{sample}-{unit}/Aligned.toTranscriptome.out.bam",
        "results/star/pe/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-pe/{sample}-{unit}.log",
    params:
        index=config["star"]["star-genome"],
        extra="--quantMode GeneCounts TranscriptomeSAM "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFilterIntronMotifs RemoveNoncanonical "
        "--chimSegmentMin 10 "
        "--chimOutType SeparateSAMold "
        "--outSAMunmapped Within "
        "--sjdbGTFfile {} {}".format(
            config["star"]["gtf"], config["star"]["params"]
        ),
    threads: 16
    resources:
        mem_mb=20000,
    envmodules:
        "star"
    wrapper:
        "v1.3.2/bio/star/align"

STAR
    --genomeDir $genomedir
    --sjdbGTFfile $GTF_mapp
    --readFilesIn $fastq1 $fastq2
    --outSAMtype BAM SortedByCoordinate
    --readFilesCommand zcat
    --alignIntronMin 30 
    --alignIntronMax 1200000
    --runThreadN 32

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
        "results/star/{ends}/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools/index/{sample}-{unit}.{ends}.sortedByCoord.log"
    resources:
        mem_mb=20000,
    envmodules:
        "samtools"
    shell:
        """
        samtools index {input[0]} 2>&1 {log}
        """
