# Snakemake workflow: rna-seq-star-deseq2

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737358.svg)](https://doi.org/10.5281/zenodo.4737358)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/rna-seq-star-deseq2/workflows/Tests/badge.svg?branch=master)](https://github.com/snakemake-workflows/rna-seq-star-deseq2/actions?query=branch%3Amaster+workflow%3ATests)


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about">About</a>
    </li>
    <li>
      <a href="#updates">Updates</a>
    </li>
    <li>
      <a href="#quick-setup">Quick Setup</a>
    </li>  
    <li>
      <a href="#dev-setup">Dev Setup</a>
    </li>
    <li>
      <a href="#workflow-dag">Workflow-DAG</a>
    </li>
    <li>
      <a href="#known-issues">Known-Issues</a>
    </li>
  </ol>
</details>




## About

This workflow performs a differential gene expression analysis with STARv2, RSEM and Deseq2; parameters and workflow establisehd by the PughLab (https://github.com/pughlab).
It has been customized to work on the HPC4Health Slurm cluster and includes extra analysis to do genotype matching.

Link to the original snakemake workflow: [snakemake-workflows/rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)

## Updates
* [2021-10-05]: (v1.0.0) Added first functional RNAseq workflow, from fastq files to STAR alignment, RSEM counts, DESeq2 diffexp, and clusterProfiler GSEA and over-representation analysis
* [TODO]: Add a script to initiate a blank working copy for building envs
* [TODO]: Implemented first working copy of genotype identification

## Quick-Setup
As states above, this workflow is set up to work on the H4H cluster. As such, there are major limitations when working on this system. For instance, only the home directory contains internet access while typically, all analysis needs to be conducted on the project directories which can only be accessed on nodes with no internet access. Additionally, the home directory contains a very limited amount of storage space.

### 0. Install Snakemake as a conda env (Build node)
More detailed information can be found on the [snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### 1. Clone the workflow (Build node)
We first want to clone the workflow to our workflow directory in our home. This requires internet access, so be sure to log in to the build node first and activate your snakemake conda env:
```
cd ~/workflows
git clone git@github.com:mcgahalab/rna-seq-star-deseq2.git
```

### 2. Setup your project directory (interactive-node)
Next we want to set up the project directory, this requires some symlinking and formatting of the file names.
Within your project directory, you will need to create the following directories:
```
mkdir config data logs resources results
```

and we want to symlink the following directories:
```
workflowpath='/path/to/rna-seq-star-deseq2/
ln -s ${workflowpath}/workflow/scripts .
ln -s ${workflowpath}/slurm .
```

#### 2.a) Setup the data directory
While the data directory shouldnt have any limitations on the sample naming convention, my files usually take on the following format:
```
[sample].R1.merged.fastq.gz  #single-end
[sample].R1.merged.fastq.gz  #paired-end
[sample].R2.merged.fastq.gz  #paired-end
e.g. DAB-C-1_S10.R1.merged.fastq.gz 
```

#### 2.b) Setup the config directory
The config directory contains 3 major files:

  * `config.yaml`: Contains all the parameters for running the workflow, specific to the project. Make sure to set up all the absolute paths to the reference files where appropriate. If you want to test multiple conditions for diffexp or alter the deseq2 modelling, change the model parameter. The conditions for testing need to be consistent with your `samples.tsv`:
```
diffexp:
  # contrasts for the deseq2 results method
  contrasts:
    dabc-vs-dabg:
      - DAB-C
      - DAB-G
    dmsoc-vs-dmso-g:
      - DMSO-C
      - DMSO-G
  model: ~condition + date
```
  * `samples.tsv`: This file contains your sample sheet with the different conditions. This does not have to be in synch with your actual fastq file IDs:
```
sample_name	condition
DAB-C-2_S11	DAB-C
DAB-C-3_S12	DAB-C
DAB-G-1_S13	DAB-G
```
  * `units.tsv`: This file contains the relative paths to your FASTQ files. The sample names have to be consistent with your `samples.tsv` file. If single-end RNA files are given, only populate the fq1 column.
```
sample_name	unit_name	fq1	fq2	sra	adapters	strandedness

DAB-C-2_S11	merged	data/DAB-C-2_S11.R1.merged.fastq.gz	data/DAB-C-2_S11.R2.merged.fastq.gz
DAB-C-3_S12	merged	data/DAB-C-3_S12.R1.merged.fastq.gz	data/DAB-C-3_S12.R2.merged.fastq.gz
DAB-G-1_S13	merged	data/DAB-G-1_S13.R1.merged.fastq.gz	data/DAB-G-1_S13.R2.merged.fastq.gz
```

### 3. Running the workflow
With the project directory created and the configuration file, samples file, and units file all setup, we are ready to run the workflow. At this point, we need to navigate back to the workflow directory (`/path/to/rna-seq-star-deseq2`). The current `scheduler.sh` is set to use my conda environments and snakemake-wrappers. However, these are subject to change with no notice, based on ongoing development. If you want to build your own environments and wrappers, refer to the **Dev-Setup** section.  

You will need to configure your `workflows/Snakefile` file and change the `workdir:` path to your project directory.
```(base) 40373Mb [quever@node47 rna-seq-star-deseq2]$ head workflow/Snakefile 
# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#container: "docker://continuumio/miniconda3"

workdir: "/path/to/projectdir"
```

Beyond this change, you should be able to execute the workflow by activating snakemake, then checking the dry-run using `snakemake -n`.  If the dry-run looks good with no errors, submit the jobs using `sbatch scheduler.sh`.  If you need to increase memory/cpu/time requirements for certain rules, you can modify these in the `slurm/cluster.json` file.



## Dev-Setup


### 1. Follow Quick-Setup steps 0-2
Follow the aforementioned steps to setup your project directory and snakemake workflow.

### 2. Setup a blank analysis directory (Build/Project node)
The basic idea of this workflow is that we set up all the conda-envs using a set output path in our home directory. However, to run the `--create-env-only` command, we need to have a directory set up with all the files needed to execute the entire workflow. In the future, I will add a script to do this automatically, but for now, it must be made manually:

```
mkdir -p ~/workflows/intialize/rnaseq-star
cd ~/workflows/intialize/rnaseq-star

mkdir config data resources
ln -s $(readlink -f /cluster/home/selghamr/workflows/rna-seq-star-deseq2/scripts) .
ln -s $(readlink -f $/cluster/home/selghamr/workflows/rna-seq-star-deseq2/slurm) .

cp /cluster/home/selghamr/workflows/rna-seq-star-deseq2/config/* config/
```
At this step, you will want to populate your `config/units.tsv` and `config/samples.tsv` with some junk nonsensical data, then you'll want to `touch` all the files into creation in the `data/` directory. For example, if you choose not to modify the `units.tsv` or `samples.tsv` (it should be fine), you can touch the following files in your data dir:

```
cd data
touch  A.1.fastq.gz A.2.fastq.gz A2.1.fastq.gz A2.2.fastq.gz B.1.fastq.gz B.2.fastq.gz
```

While not listed here, many of your reference files (e.g. `genoma.fa`) will be in your Project directories, which you cannot access from the build node.  Be sure to reconfigure the `config/config.yaml` to point to reference files that are in the home directory (`resources/` directory). As nothing will actually be run, you can simply `touch` these files for the sake of satisfying the snakemake file checker.


### 3. Download the snakemake-wrappers (Build node)
Again, while in the build node, you can only access the home directory. Because of this, you will need to have a local copy of `snakemake-wrappers` in your home directory in order to build your `conda-envs`

This workflow requires: `v0.75.0` and `0.59.2`.

Visit the the git repo for [snakemake-wrappers](https://github.com/snakemake/snakemake-wrappers). Switch to the **tag** that you want to download (e.g. v0.75.0), then download the zip file (e.g. `wget 'https://github.com/snakemake/snakemake-wrappers/archive/refs/tags/v0.75.0.zip'`) to your snakemake-wrappers directory. Once unzipped, be sure to relabel it as `v0.75.0` or whatever tag you had it.

Basically, you will need to add a path that will have the following structure:
```
/path/to/snakemake-wrappers/[TAG]/bio/
```

### 4. Building the conda-envs (Build node)
Now that you have your "intialize" and your "project" directory set up, you can start building your conda environments required for the workflow.

Before you build your conda-envs, you'll need to switch the `workdir:` path in your `workflow/Snakefile`:
```
> workflow/Snakefile
workdir: "/path/to/workflows/intialize/rnaseq-star"
```

Once the Snakefile workdir path has been changed, you can run the entire snakemake workflow using `conda-create-envs-only`:

```
cd ~/workflows/rna-seq-star-deseq2
condaprefix=$(readlink -f .)"/.snakemake/conda"

snakemake \
--use-conda \
--use-singularity \
--conda-create-envs-only \
--conda-frontend conda \
--conda-prefix ${condaprefix} \
--wrapper-prefix 'file:////cluster/home/selghamr/snakemake_wrappers/' \
--cores 4
```
**Note:** the `--wrapper-prefix` is labelled as `file:///` with three `/`'s.  It WILL throw a warning/error during building/running the workflow, but it MUST be set this way or it will not work.

The idea behind this is that snakemake will install the conda envs to `.snakemake/conda`. It creates a hash to label that environment it builds. However, the hash is generated based on the `conda-prefix` and the hash of the `env.yaml` that it is building. As long as the `conda-prefix` and `env.yaml` remain unchanged, it will use the pre-existing environment found in `.snakemake/conda/[HASH]`

### 5. Run your workflow (Build/Project node)
Once the workflow has been set up and the environments have been created, you can finally run your pre-configured workflow on your project directory. You will need to reconfigure some paths in your `scheduler.sh` script. Be sure to activate your `snakemake` env before queuing the scheduler.

Before you run your workflow on your project directory, you'll need to switch the `workdir:` path in your `workflow/Snakefile`:
```
> workflow/Snakefile
workdir: "/path/to/group_directory/project"
```

The `scheduler.sh` script runs the following command in a 5-day long job, with its main purpose to track the job queue and job submission on the slurm cluster.
```
cd /cluster/home/selghamr/workflows/rna-seq-star-deseq2
condaprefix='/cluster/home/selghamr/workflows/rna-seq-star-deseq2/.snakemake/conda'

snakemake \
--jobs 6 \
--profile slurm \
--cluster-config slurm/cluster.json \
--conda-frontend conda \
--use-conda \
--use-singularity \
--conda-prefix ${condaprefix} \
--wrapper-prefix 'file:////cluster/home/selghamr/snakemake_wrappers/v0.77.0' \
--rerun-incomplete
```

You can submit the scheduler as followed:
```
conda activate snakemake
sbatch scheduler.sh
```

## Workflow-DAG
## Known-Issues
1. rsem-calculate-expressions keeps using single-end mode and throws errors when I use paired-end
  * The snakemake-wrapper for `calculate-expression` looks for a param called `paired-end` in your rule. However, the "-" makes this an incompatible param name and thus is unusable.  Instead, you need to change the `wrapper.py` for the calculate-expressions wrapper to look for `paired_end` instead. Thus, change `line 19` in `/path/to/snakemake-wrappers/0.77.0/bio/rsem/calculate-expression/wrapper.py` from paired-end to paired_end.  The rule set up right now uses the paired_end param already.
  
2. rsem-calculate-expressions throws an rsem-tbam2gbam error about two mates aligning to two different transcripts after ~1-2 hrs of runtime:
  * I don't know the source of this error, but it seems to be specific to rsem=3.3.3 from conda and fixed in the 3.3.0 version (https://github.com/deweylab/RSEM/issues/126).  The rsem-calculate-expressions environment.yaml uses the rsem=3.3.3 and will require manually editting the file to change the rsem file `/path/to/snakemake-wrappers/0.77.0/bio/rsem/calculate-expression/environment.yaml`. However, to date, I haven't been able to give rsem=3.3.0 installed in my conda env due to a glibc=2.17.0 error.  As such, there is a temporary fix set up where I use the rsem=3.3.0 installed on H4H, accessed via module load, to run that rule.

3. CreateCondaEnvironmentException related to multiqc/environment.yaml:
  * This issue is related to the bioconda installation of `mutliqc==1.7`. The simplest method to fix this issue is to change the environment.yaml of `/path/to/snakemake-wrappers/0.77.0/bio/multiqc/environment.yaml` to place conda-forge ahead of bioconda. Alternatively, you can modify the associated yaml file in `/path/to/rna-seq-star-deseq2/.snakemake/conda/XYZ.yaml` to place conda-forge ahead of the bioconda env, and then force the installation using the following code, and then change the yaml file back to its original state:
```
id='XYZ'   #e.g. '1943aee1fff2ccaea034ffd15210faca'
targetdir=$(readlink -f .)

target_env_file=${targetdir}'/'${id}'.yaml'
env_path=${targetdir}'/'${id}
mkdir ${env_path}
conda env create --file ${target_env_file} --prefix ${env_path}
touch ${env_path}/env_setup_done ${env_path}/env_setup_start
```
