#rule prepare_reference:
#  input:
#    reference_genome=config['ref_index']['genome'],
#  output:
#    seq="ref/reference.seq",
#    grp="ref/reference.grp",
#    ti="ref/reference.ti",
#  params:
#    extra="--gtf {}".format(config["star"]["gtf"]),
#  log:
#    "logs/rsem/prepare-reference.log",
#  wrapper:
#    "0.77.0/bio/rsem/prepare-reference"

rule prepare_reference:
  input:
    reference_genome=config['ref_index']['genome']
  output:
    seq=config['ref_index']['rsem'] + "/reference.seq",
    grp=config['ref_index']['rsem'] + "/reference.grp",
    ti=config['ref_index']['rsem'] + "/reference.ti"
  params:
    refpath=config['ref_index']['rsem'],
    # extra="--gtf {}".format(config["star"]["gtf"]),
    extra="--gtf " + base_ref + "/" + ref_annots,
    scripts_path=BASE_SCRIPTS
  log:
    "logs/rsem/prepare-reference.log"
  threads: 1
  resources:
    mem_mb=2000
  envmodules:
    "rsem"
  shell:
    """
    # if [ -f "{params.refpath}.seq" ]; then
    #   ln -s {params.refpath}.* ref/
    # else
    #module load rsem/1.3.0
    python {params.scripts_path}/prepare-rsem-reference.py {input.reference_genome} {output.seq} '{params.extra}' {log}
    # fi
    """

#rule calculate_expression:
#  input:
#    bam=get_star_transcriptome,
#    reference="ref/reference.seq",
#  output:
#    genes_results="results/rsem/{sample}-{unit}.genes.results",
#    isoforms_results="results/rsem/{sample}-{unit}.isoforms.results",
#    tbam=temp("results/rsem/{sample}-{unit}.transcript.bam"),
#    gbam=temp("results/rsem/{sample}-{unit}.genome.bam"),
#  params:
#    paired_end=lambda w: is_paired_end(w.sample),
#    extra="--estimate-rspd --output-genome-bam --time --forward-prob 0 --seed 42",
#  log:
#    "logs/rsem/calculate_expression/{sample}-{unit}.log",
#  wrapper:
#    "0.77.0/bio/rsem/calculate-expression"

rule calculate_expression:
  input:
    bam=get_star_transcriptome,
    # reference="ref/reference.seq",
    reference=rules.prepare_reference.output[0],
  output:
    genes_results=BASE_OUT+ "/2.RSEM/{sample}-{unit}.genes.results",
    isoforms_results=BASE_OUT+ "/2.RSEM/{sample}-{unit}.isoforms.results",
    tbam=temp(BASE_OUT+ "/2.RSEM/{sample}-{unit}.transcript.bam"),
    gbam=temp(BASE_OUT+ "/2.RSEM/{sample}-{unit}.genome.bam"),
  params:
    outprefix=BASE_OUT + "/2.RSEM/{sample}-{unit}",
    paired_end=lambda w: "--paired-end" if is_paired_end(w.sample) else "",
    extra="--bam --estimate-rspd --output-genome-bam --time --forward-prob 0 --seed 42",
  threads: 2
  resources:
    mem_mb=5000
  log:
    "logs/rsem/calculate_expression/{sample}-{unit}.log",
  envmodules:
    "rsem"
  shell:
    """
    ref=$(echo {input.reference} | sed 's/\\..*//');
    # module load rsem/1.3.0;
    rsem-calculate-expression --num-threads {threads} {params.extra} {params.paired_end} {input.bam} $ref {params.outprefix} > {log} 2>&1
    """

#rule rsem_generate_data_matrix:
#  input:
#    get_rsem_output_all_units,
#  output:
#    temp("results/counts/all.tmp"),
#  params:
#    extra="",
#  log:
#    "logs/rsem/generate_data_matrix.log",
#  wrapper:
#    "0.77.0/bio/rsem/generate-data-matrix"

rule rsem_generate_data_matrix:
  input:
    get_rsem_output_all_units,
  output:
    temp(BASE_OUT+ "/2.RSEM/counts/all.tmp"),
    # temp("results/counts/all.tmp"),
  params:
    extra="",
  log:
    "logs/rsem/generate_data_matrix.log",
  threads: 1
  resources:
    mem_mb=5000
  envmodules:
    "rsem"
  shell:
    """
    # module load rsem/1.3.0;
    rsem-generate-data-matrix {params.extra} {input} > {output} 2>{log}
    """

## Need to make this rule more robust for non-merged unit0-names
rule format_data_matrix:
  input:
    BASE_OUT+ "/2.RSEM/counts/all.tmp"
    # "results/counts/all.tmp",
  output:
    # "results/counts/all.tsv",
    BASE_OUT+ "/2.RSEM/counts/all.tsv",
  log:
    "logs/rsem/format_data_matrix.log",
  threads: 1
  resources:
    mem_mb=5000
  shell:
    """
    sed 's/\"//g' {input} | sed 's/results\/rsem\///g' | sed 's/-merged.genes.results//g' | sed '1 s/^/gene/' | sed -e 's/\.[0-9]*//g' > {output}
    """

rule rsem_generate_tpm_matrix:
  input:
    get_rsem_output_all_units,
  output:
    # temp("results/counts/all_tpm.tmp"),
    temp(BASE_OUT+ "/2.RSEM/counts/all_tpm.tmp"),
  params:
    extra="TPM",
  log:
    "logs/rsem/rsem_generate_tpm_matrix.log",
  threads: 1
  resources:
    mem_mb=5000
  envmodules:
    "rsem"
  shell:
    """
    # module load rsem/1.3.0
    perl scripts/rsem-generate-data-matrix-modified.pl {params.extra} {input} > {output} 2>{log}
    """

rule format_tpm_matrix:
  input:
    # "results/counts/all_tpm.tmp",
    BASE_OUT+ "/2.RSEM/counts/all_tpm.tmp"
  output:
    # "results/counts/all_tpm.tsv",
    BASE_OUT+ "/2.RSEM/counts/all_tpm.tsv"
  log:
    "logs/rsem/format_tpm_matrix.log",
  threads: 1
  resources:
    mem_mb=5000    
  shell:
    """
    sed 's/\"//g' {input} | sed 's/results\/rsem\///g' | sed 's/-merged.genes.results//g' | sed '1 s/^/gene/' | sed -e 's/\.[0-9]*//g' > {output}
    """
