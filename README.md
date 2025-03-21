# Hevesi_2023

1. System requirements
 - All software dependencies and operating systems (including version numbers) are available in the [methods section][methods]
 - Archived version of this repo has been tested on x86_64-pc-linux-gnu (64-bit) platform running under Ubuntu 22.04.1 LTS with R version 4.2.2 and Python version 3.8.8
 - GPU is required to run ambient RNA removal with CellBender
2. Installation guide
 - Installation is available via Snakemake framework infrastructure. The Snakefile pipeline is in the project root. Before run: 1) You need to edit path to the cellranger environment file on your computer ($CELLRANGER). 2) Be sure that Apptainer (former Singularity is installed and available in $PATH).  
 - Typical install time on a "normal" desktop computer should be within a hour that you need to install cellranger, Anaconda, Apptainer/Singularity, to create snakemake conda environment and to download all configured container-images.
3. Demo
 - We include preprocessed data that one can use to reproduce panels of the figure 2 related to snRNA-seq data from the manuscript.
 - It is expected that output images can slightly differ due to random factor native to dimension reduction technics.
 - Expected run time for demo on a "normal" desktop computer is within a hour.
4. Instructions for use
 - To run the demo on preprocessed data you need to source code/render.R file, which will rerun rendering of Rmarkdown analysis notebooks for figure 2 via workflowr package.
 - To reproduce whole pipeline you would need to modify and run command `snakemake -f --use-singularity -prk --cores <FILL-IN-N-CORES> --resources nvidia_gpu=1 --singularity-args "--nv --bind <FILL-IN-PATH-TO-CLONED-REPO>" --rerun-incomplete --retries 1 --resources mem_mb=<FILL-IN-MEM>`

[methods]: https://harkany-lab.github.io/Hevesi_2023/methods.html "Methods section of the analysis website"
