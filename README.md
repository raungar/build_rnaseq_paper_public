# Pipeline to implement quantification and outlier detection across genome builds
**Contact:** Stephen Montgomery (smontgom@stanford.edu), Rachel Ungar (raungar@stanford.edu), Page Goddard (pgoddard@stanford.edu), Tanner Jensen (tannerj@stanford.edu)
This pipeline is an iteration on the pipeline released in the paper [Identification of rare-disease genes using blood transcriptome sequencing and large control cohorts](https://www.nature.com/articles/s41591-019-0457-8)



# Pipeline Overview
![alt text](GenomeBuildPipeline.png)

# Installation Preparation
* if necessary, install Miniconda as recommended [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

    ```
    wget https://docs.conda.io/en/latest/miniconda.html#linux-installers
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

* install snakemake via conda and mamba as recommended [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

    ```
    conda install -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -n snakemake snakemake
    conda activate snakemake
    snakemake --help
    ```

* Set up a snakemake slurm profile (good install walkthrough [here](http://bluegenes.github.io/Using-Snakemake_Profiles/)) and information on actual use [here](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/)
    * Setting up the slurm profile requires cookiecutter, which can be installed with conda:

    ```
    conda install -c conda-forge cookiecutter
    ```

# environment_paths.yaml
There should be no direct paths in the snakemake. Please add any paths to the yaml file.