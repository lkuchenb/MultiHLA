# HLA Typing Workflow

<p align="center">
<img src="img/multihla_diag.svg?raw=true" alt="Workflow Diagram"/>
</p>

## Scope of this workflow

This workflow enables the concurrent analysis of WES or WGS data using
publicly available software to derive HLA haplotypes from this type of data.

## Currently available software tools

 * xHLA

	*Xie, C., Yeo, Z. X., Wong, M., Piper, J., Long, T., Kirkness, E. F., ... & Brady, C. (2017). Fast and accurate HLA typing from short-read next-generation sequence data with xHLA. Proceedings of the National Academy of Sciences, 114(30), 8059-8064.*

	The workflow implements read mapping the reads against hg38 without alt
	contigs using `bwa mem` as instructed by the authors. The mapped reads are then
	sorted and index using samtools.

	The workflow utilizes the [Docker Image][2] provided by the authors to
	perform the actual HLA typing.

 * HLA-VBSeq

	*Nariai, N., Kojima, K., Saito, S., Mimori, T., Sato, Y., Kawai, Y., ... & Nagasaki, M. (2015, December). HLA-VBSeq: accurate HLA typing at full resolution from whole-genome sequencing data. In BMC genomics (Vol. 16, No. S2, p. S7). BioMed Central.*

	*Wang, Y. Y., Mimori, T., Khor, S. S., Gervais, O., Kawai, Y., Hitomi, Y., ... & Nagasaki, M. (2019). HLA-VBSeq v2: improved HLA calling accuracy with full-length Japanese class-I panel. Human Genome Variation, 6(1), 1-5.*

	The workflow implements read mapping the reads against hg19 without alt
	contigs. The authors instructions merely state to "map against hg19"
	without any further specifics, but mapping against hg19 with alt contigs
	yielded very poor typing results with missing HLA class I genes, thus
	the workflow uses hg19 without alt contigs.

	HLA-VBSeq released two reference database versions:

	* v1 database based on IMGT/HLA database, Release 3.15.0
	* v2 database based on IMGT/HLA database Release 3.31.0 and Japanese HLA reference dataset

 * OptiType

	*Szolek, A., Schubert, B., Mohr, C., Sturm, M., Feldhahn, M., & Kohlbacher, O. (2014). OptiType: precision HLA typing from next-generation sequencing data. Bioinformatics, 30(23), 3310-3316.*

	The workflow invokes the [OptiType snakemake wrapper][1] without prior filtering
	of reads.

* HLA-LA

	*Dilthey, A. T., Mentzer, A. J., Carapito, R., Cutland, C., Cereb, N., Madhi, S. A., ... & Phillippy, A. M. (2019). HLA* LA—HLA typing from linearly projected graph alignments. Bioinformatics, 35(21), 4394-4396.*

	The workflow uses reads mapped against the human genome (hg38) without
	alt contigs as input for HLA-LA. A corresponding reference txt file for HLA-LA
	is part of this workflow repository. The preprocessed graph directory
	`PRG_MHC_GRCh38_withIMGT` can be either placed manually in
	`typing/hla_la/hla_la.graphs/` or it will be downloaded and preprocessed
	automatically.

	The workflow uses the [HLA-LA bioconda package][3] for graph preprocessing and HLA typing.


## Usage

 1. Install snakemake

	```
	conda install -c conda-forge mamba
	mamba create -c conda-forge -c bioconda -n snakemake snakemake
	conda activate snakemake
	```

 1. Clone the *MultiHLA* repository
	```
	git clone https://github.com/lkuchenb/MultiHLA.git hla_typing
	cd hla_typing
	```

 1. Put the input files in place<br/>
    *MultiHLA* comes with a predefined folder structure:
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

		Place or link the required human genome references here as described for each supported method, otherwise they will be automatically downloaded.
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
	process all samples of a dataset, for example the dataset `dataset_1`
	described in `datasets/dataset_1.tsv` use
	```
	snakemake --use-conda --use-singularity typing/dataset_1.all.multihla
	```
	Memory and run time requirements for each job are noted in their resources (`mem_mb` and `time`).

[1]: https://snakemake-wrappers.readthedocs.io/en/0.61.0/wrappers/optitype.html
[2]: https://hub.docker.com/r/humanlongevity/hla
[3]: https://bioconda.github.io/recipes/hla-la/README.html
