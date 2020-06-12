# HLA Typing Benchmark Workflow

<p align="center">
<img src="img/hla_bench_diag.svg?raw=true" alt="Workflow Diagram"/>
</p>

## Scope of this workflow

This workflow enables the concurrent analysis of WES or WGS data using
publicly available software to derive HLA haplotypes from this type of data.

## Currently available software tools

 * xHLA

	The workflow implements read mapping the reads against hg38 without alt
	contigs using `bwa mem` as instructed by the authors. The mapped reads are then
	sorted and index using samtools and the HLA typing is performed using the
	docker container provided by the authors.

	The workflow requires the human genome reference `hg38` with no alt
	contigs and with index files produced by `bwa index` to be provided in the
	`ref/` folder under the following names:
	```
	hg38.noalt.fa
	hg38.noalt.fa.amb
	hg38.noalt.fa.ann
	hg38.noalt.fa.bwt
	hg38.noalt.fa.pac
	hg38.noalt.fa.sa
	```

## Usage

 1. Install snakemake

	```
	conda install -c conda-forge mamba
	mamba create -c conda-forge -c bioconda -n snakemake snakemake
	conda activate snakemake
	```

 1. Clone the *HLABench* repository
	```
	git clone https://github.com/lkuchenb/HLABench.git hla_benchmark
	cd hla_benchmark
	```

 1. Put the input files in place<br/>
    *HLABench* comes with a predefined folder structure:
    * `dataset/`

		A dataset is defined as a set of samples. Place a TSV file here for every dataset with the following three named columns:
		```
		SampleName  FileNameR1                              FileNameR2
		Donor1      SEQ_D1_DAT_01_S53_L001_R1_001.fastq.gz  SEQ_D1_DAT_01_S53_L001_R2_001.fastq.gz
		Donor1      SEQ_D1_DAT_01_S53_L002_R1_001.fastq.gz  SEQ_D1_DAT_01_S53_L002_R2_001.fastq.gz
		Donor2      SEQ_D2_DAT_01_S54_L001_R1_001.fastq.gz  SEQ_D2_DAT_01_S54_L001_R2_001.fastq.gz
		Donor2      SEQ_D2_DAT_01_S54_L002_R1_001.fastq.gz  SEQ_D2_DAT_01_S54_L002_R2_001.fastq.gz
		Donor3      SEQ_D3_DAT_01_S55_L001_R1_001.fastq.gz  SEQ_D3_DAT_01_S55_L001_R2_001.fastq.gz
		Donor3      SEQ_D3_DAT_01_S55_L002_R1_001.fastq.gz  SEQ_D3_DAT_01_S55_L002_R2_001.fastq.gz
		```
		FASTQ files have to come in gziped pairs and be named
		`{prefix}_R[12]{suffix}.fastq.gz`. A sample can be covered by an arbitrary
		number of FASTQ pairs (at least one).
    * `fastq/`

		Place the FASTQ files as listed in your dataset sheet here.
    * `ref/`

		Place the required human genome references here as described for each supported method.
    * `trim/`

		This is an output folder. It will be filled with adapter trimmed versions of the provided FASTQ files.
    * `typing/{method}/`

		This is an output folder. It will be filled with subfolders for each method.
    * `workflow/`

		This folder contains the workflow code.
 1. Run the workflow

	Invoke snakemake using `snakemake --use-conda --use-singularity`. This enables
	snakemake to automatically install dependencies into conda environments that
	are created on the fly and also enables the container based jobs to run. To
	process all samples of a dataset, for example the dataset `benchmark_1`
	described in `datasets/benchmark_1.tsv` use
	```
	snakemake --use-conda --use-singularity typing/xhla/benchmark_1.tsv
	```
	Memory and run time requirements for each job are noted in their params (`cluster_mem` and `cluster_rt`). For example, to interface with an SGE grid engine using DRMAA you may use
	```
	snakemake --latency-wait 20 --drmaa " -V -l h_vmem={params.cluster_mem} -l h_rt={params.cluster_rt} -pe smp {threads} -j yes -o ${PWD}/cluster_log" -j 999 --use-conda --use-singularity
	```
