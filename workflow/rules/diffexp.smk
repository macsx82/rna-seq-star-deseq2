rule count_matrix:
    input:
        get_star_output_all_units,
    output:
        # "results/counts/all_old.tsv",
        BASE_OUT+"/3.DESeq2/all_old.tsv"
    log:
        "logs/count-matrix.log",
    threads: 1
    resources:
        mem_mb=5000
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
        scripts_path=BASE_SCRIPTS
    # conda:
    #     "../envs/pandas.yaml"
    envmodules:
        "python/3.10.3"
    script:
        "{params.scripts_path}/count-matrix.py"


rule deseq2_init:
    input:
        # counts="results/counts/all.tsv",
        #set it to work with DESeq2 COUNTS instead of RSEM
        counts=rules.count_matrix.output[0]
        #point to the usage of the RSEM matrix data
        # counts=rules.format_data_matrix.output[0]
    output:
        BASE_OUT+"/3.DESeq2/all.rds",
        BASE_OUT+"/3.DESeq2/normcounts.tsv",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
        scripts_path=BASE_SCRIPTS
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    resources:
        mem_mb=20000
    script:
        "{params.scripts_path}/scripts/deseq2-init.R"


# rule pca:
#     input:
#         "results/deseq2/all.rds",
#     output:
#         report("results/pca.pdf", "../report/pca.rst"),
#     params:
#         pca_labels=config["pca"]["labels"],
#     conda:
#         "../envs/deseq2.yaml"
#     log:
#         "logs/pca.log",
#     script:
#         "../scripts/plot-pca.R"


# rule deseq2:
#     input:
#         "results/deseq2/all.rds",
#     output:
#         table=report(
#             "results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"
#         ),
#         ma_plot=report("results/diffexp/{contrast}.ma-plot.pdf", "../report/ma.rst"),
#     params:
#         contrast=get_contrast,
#         species=config['ref']['species'],
#     conda:
#         "../envs/deseq2.yaml"
#     log:
#         "logs/deseq2/{contrast}.diffexp.log",
#     threads: get_deseq2_threads
#     script:
#         "../scripts/deseq2.R"

# rule gsea:
#     input:
#         "results/diffexp/{contrast}.diffexp.tsv",
#     output:
#         ora_tbl=report("results/gsea/{contrast}.ora.tsv", "../report/ora-go.rst"),
#         gsea_tbl=report("results/gsea/{contrast}.gsea.tsv", "../report/gsea-kegg.rst"),
#         gsea_pdf="results/gsea/{contrast}.gsea.pdf",
#     params:
#         genome=config['ref']['species'],
#         contrast=get_contrast,
#         minbase=config['gsea']["min_base_mean"],
#         maxp=config['gsea']["max_p"],
#     conda:
#         "../envs/deseq2.yaml"
#     log:
#         "logs/gsea/{contrast}.log",
#     threads: get_deseq2_threads
#     script:
#         "../scripts/gsea.R"
