# PopGen CLI Iterative Gvcf Genotyper User Documentation

## Table of Contents

1. [Installation](#installation)
   1. [Get ICA API key](#get-ica-api-key)
   2. [Download and install ICA CLI](#download-and-install-ica-cli)
   3. [Create demo project and link the PopGen CLI bundle](#create-demo-project-and-link-the-popgen-cli-bundle)
   4. [Download and install PopGen CLI](#download-and-install-popgen-cli)
   5. [Install the PopGen CLI in a virtual environment](#install-the-popgen-cli-in-a-virtual-environment)
2. [Launch the IGG demo](#launch-the-igg-demo)
   1. [Download IGG demo](#download-igg-demo)
   2. [Configure IGG demo](#configure-igg-demo)
   3. [Upload the configuration](#upload-the-configuration)
   4. [Submit demo analyses](#submit-demo-analyses)
      1. [Submit version-1 IGG analyses](#submit-version-1-igg-analyses)
      2. [Submit version-2 IGG analyses](#submit-version-2-igg-analyses)
      3. [Monitor IGG analyses](#monitor-igg-analyses)
3. [IGG output ICA projects and data](#igg-output-ica-projects-and-data)
   1. [Output data of step1](#output-data-of-step1)
   2. [Version specific output data](#version-specific-output-data)
   3. [Output data of step2](#output-data-of-step2)
   4. [Output data of step3](#output-data-of-step3)
   5. [Output data of step4](#output-data-of-step4)
   5. [Output data of step5](#output-data-of-step5)
4. [Details about the IGG demo](#details-about-the-igg-demo)
   1. [Batch](#batch)
   2. [Version](#version)
   3. [Gender information](#gender-information)
   4. [gVCF version](#gvcf-version)
   5. [Shard](#shard)
   6. [Subshard](#subshard)
   7. [Chromosome](#chromosome)
   8. [Reference Genome](#reference-genome)
   9. [Annotation](#annotation)
   10. [Concatenation by output format](#concatenation-by-output-format)

## Installation

### Get ICA API key

Follow the instructions in [Generate API key](https://help.ica.illumina.com/get-started/gs-getstarted#generate-an-api-key) to get an API key.

### Download and install ICA CLI

Follow the instructions in [Install ICA CLI](https://help.ica.illumina.com/command-line-interface/cli-installation) to install the ICA CLI.

After installation, follow the instructions in [Configure CLI](https://help.ica.illumina.com/command-line-interface/cli-authentication) to configure the CLI with the API key.

```bash
icav2 config set
```

### Create demo project and link the PopGen CLI bundle

PopGen CLI is currently available in the following ICA regions.

| region | region code | country |
|--------|-------------|---------|
| use1   | US          | United States |
| usw2   | US-OR       | United States |
| cac1   | CA          | Canada      |
| euw2   | UK          | United Kingdom |
| euc1   | EU          | Germany     |
| apne1  | JP          | Japan       |
| apne2  | KR          | Republic of Korea |
| apse1  | SG          | Singapore   |
| apse2  | AU          | Australia   |
| apse3  | ID          | Indonesia   |
| aps1   | IN          | India       |
| ilc1   | IL          | Israel      |
| mec1   | AE          | United Arab Emirates |

Check if you already have a demo project called `<region>-popgen-cli-release-demo` in your domain. You may check this either within the [Web UI](https://ica.illumina.com/ica/projects) or use the command line

```bash
icav2 projects list
```

If the project does not exist, create one with the name  `<region>-popgen-cli-release-demo`. For details, see [Create new project](https://help.ica.illumina.com/home/h-projects#create-new-project).

Link the bundle `PopGen CLI (<region code>)` to the demo project. For details of linking the bundle (https://help.ica.illumina.com/home/h-bundles#linking-an-existing-bundle-to-a-project)

### Download and install PopGen CLI

Once the bundle is linked to the release demo project, you may find the PopGen CLI release under folder `/popgen-cli-release/v<version>/`. The latest version is `2.2.1`.

Download the PopGen CLI python wheel file from Web UI or using ICA CLI

```bash
icav2 projects enter <region>-popgen-cli-release-demo
icav2 projectdata download /popgen-cli-release/v<version>/popgen_cli-<version>-py3-none-any.whl
```

### Install the PopGen CLI in a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate
pip3 install popgen_cli-<version>-py3-none-any.whl
popgen-cli --version
```

## Launch the IGG demo

To launch IGG demo, there are 4 steps:

1. Download IGG demo
2. Configure IGG demo
3. Upload the configuration
4. Submit demo analyses

### Download IGG demo

Download the following files related to the IGG demo

- the demo tar file `popgen-dragen-igg-demo.tar.gz`
- the README file about IGG usage and demo `README.dragen-igg.md`

```bash
icav2 projectdata download /popgen-cli-release/v<version>/popgen-dragen-igg-demo.tar.gz
icav2 projectdata download /popgen-cli-release/v<version>/README.dragen-igg.md
tar xvzf popgen-dragen-igg-demo.tar.gz
cd popgen-dragen-igg-demo
```

The demo folder `popgen-dragen-igg-demo` contains following scripts and folder

| Demo file/folder | Description |
|-----------------|-------------|
| 1.config.sh | the script to configure IGG project |
| 2.upload.sh | the script to upload local demo project to ICA project |
| 3.submit.sh | the script to launch IGG analyses in ICA |
| poll.sh | the script to poll ICA analyses status |
| project-data | the folder with demo metadata and configurations |
| README.dragen-igg.md | this user guide in markdown format |

### Configure IGG demo

Run following script to configure the IGG project.

```bash
./1.config.sh
```

This configuration step collects static configuration metadata that is required to run every IGG analysis. It also creates a job project in ICA, where IGG metadata is stored, pipelines are deployed, analyses are managed, and analysis log files are stored.

For example, the job project name may look like `<region>-dragen-igg-<timestamp>-jobs`. IGG output data will be stored in related projects

- `<region>-dragen-igg-<timestamp>-cohort-census`: this project stores the cohort and census files in all batches and shards (output from step1). It will be created by IGG `step1` analyses.
- `<region>-dragen-igg-<timestamp>-msvcf-version-<version-id>`: this project stores, for each aggregation version, the output msVCF, PGEN and annotation data (output from `step2`, `step3`, `step4` and `step5`). This project will be created by IGG `step2` analyses.

The result of the configuration step is stored locally in the file `project-data/secret/project-config.<timestamp>.json`.

### Upload the configuration

Run the following script to upload all metadata under `project-data` to ICA IGG job project `<region>-dragen-igg-<timestamp>-jobs` under `/meta/` folder.

```bash
./2.upload.sh
```

### Submit demo analyses

The IGG demo consists of two versions of aggregation of a total of 7 samples on chromosome 20 (`shard-88` and `shard-89`).

- Version 1: 4 samples (`batch-1`)
- Version 2: 4 samples (`batch-1`) + 3 samples (`batch-2`)

#### Submit `version-1` IGG analyses

The analyses of each step need to be completed (in `SUCCEEDED` status) before submitting the analyses in the next step.

```bash
batch_ids=1
shard_ids=88-89
version_id=1
chrom_ids=20

# step1 aggregate gvcf (batch-1, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-cohort-census
./3.submit.sh step1 $batch_ids $shard_ids

# step2 aggregate census, run annotation (version-1, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-1
./3.submit.sh step2 $version_id $shard_ids

# collect subshard information (version-1, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-1
./3.submit.sh count_subshards $version_id $shard_ids

# step3 generate msVCF, ML filtering, PGEN (version-1, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-1
./3.submit.sh step3 $version_id $shard_ids

# step4 concat msVCF, PGEN (version-1, chrom-20)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-1
./3.submit.sh step4 $version_id $chrom_ids

# step5 annotate site VCF (version-1, chrom-20)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-1
./3.submit.sh step5 $version_id $chrom_ids
```

####Submit version-2 IGG analyses

After `version-1` IGG analyses are completed, let us submit analyses for `version-2` to iteratively aggregate 3 more samples in `batch-2` in addition to the 4 samples in `batch-1`.

For `step1`, only `batch-2` analyses are required to submit. The output of `version-2` analyses (`step2`, `step3`, `step4`, `step5`) will be stored in a separate project, different from the `version-1` output project.

The analyses of each step need to be completed (in `SUCCEEDED` status) before submitting the analyses in the next step.

```bash
batch_ids=2
shard_ids=88-89
version_id=2
chrom_ids=20

# step1 aggregate gvcf (batch-2, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-cohort-census
# NOTE: no need to resubmit analyses on batch-1, shard-88 and shard-89
./3.submit.sh step1 $batch_ids $shard_ids

# step2 aggregate census, run annotation (version-2, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-2
./3.submit.sh step2 $version_id $shard_ids

# collect subshard information (version-2, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-2
./3.submit.sh count_subshards $version_id $shard_ids

# step3 generate msVCF, ML filtering, PGEN (version-2, shard-88 and shard-89)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-2
./3.submit.sh step3 $version_id $shard_ids

# step4 concat msVCF, PGEN (version-2, chrom-20)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-2
./3.submit.sh step4 $version_id $chrom_ids

# step5 annotate site VCF (version-2, chrom-20)
# output project: <region>-dragen-igg-<timestamp>-msvcf-version-2
./3.submit.sh step5 $version_id $chrom_ids
```

#### Monitor IGG analyses

You may monitor the status of analyses from ICA Web UI of the IGG job project `<region>-dragen-igg-<timestamp>-jobs`, or simply run the script

```bash
./poll.sh
```

## IGG output ICA projects and data

### Output data of `step1`

The output data of `step1` is stored in the project `<region>-dragen-igg-<timestamp>-cohort-census`. The files are split by shard and batch. The main output files, `cohort` and `census` files, are located in the following folders.

```bash
/data/batch-1/shard-88/dragen.cht.gz
/data/batch-1/shard-88/dragen.cns.gz
/data/batch-1/shard-89/dragen.cht.gz
/data/batch-1/shard-89/dragen.cns.gz
/data/batch-2/shard-88/dragen.cht.gz
/data/batch-2/shard-88/dragen.cns.gz
/data/batch-2/shard-89/dragen.cht.gz
/data/batch-2/shard-89/dragen.cns.gz
```

### Version specific output data

The output data of `step2`, `step3`, `step4`, and `step5` is per version and is stored in the project `<region>-dragen-igg-<timestamp>-msvcf-version-<version>`.

### Output data of `step2`

The output data of `step2` is stored in the project `<region>-dragen-igg-<timestamp>-msvcf-version-<version>` under the folder `/data/global-census`. The files are stored per shard.

The main output files, `global census`, multi-allelic sites VCFs are located in the following folders.

```bash
/data/global-census/shard-88/dragen.cns.gz
/data/global-census/shard-88/dragen.sites.vcf.gz
/data/global-census/shard-89/dragen.cns.gz
/data/global-census/shard-89/dragen.sites.vcf.gz
```

The number of subshards per shard and the subshard region definition and statistics are stored in the following files.

```bash
/data/global-census/shard-88/dragen.shard.num-subshards.csv
/data/global-census/shard-88/dragen.subshards.stats.tsv
/data/global-census/shard-89/dragen.shard.num-subshards.csv
/data/global-census/shard-89/dragen.subshards.stats.tsv
```

### Output data of `step3`

The output data of `step3` is stored in the project `<region>-dragen-igg-<timestamp>-msvcf-version-<version>` under the folder `/data/shard-msvcf`. The files are located per shard and subshard.

For each subshard, the full multi-sample VCF file (msVCF) is located in the following folders.

```bash
/data/shard-msvcf/shard-88/subshard-1/dragen.vcf.gz
/data/shard-msvcf/shard-89/subshard-1/dragen.vcf.gz
```

The post-processed reduced msVCF files (biallelic with cohort-level ML filtering), the site VCF files, and subshard statistics files are located in the following folders.

```bash
/data/shard-msvcf/shard-88/subshard-1/postproc/vcf/dragen.vcf.gz
/data/shard-msvcf/shard-88/subshard-1/postproc/vcf/dragen.sites.vcf.gz
/data/shard-msvcf/shard-88/subshard-1/postproc/vcf/dragen.subshards.stats.tsv
/data/shard-msvcf/shard-89/subshard-1/postproc/vcf/dragen.vcf.gz
/data/shard-msvcf/shard-89/subshard-1/postproc/vcf/dragen.sites.vcf.gz
/data/shard-msvcf/shard-89/subshard-1/postproc/vcf/dragen.subshards.stats.tsv
```

The post-processed PGEN files are located in the following folders.

```bash
/data/shard-msvcf/shard-88/subshard-1/postproc/pgen/dragen.pgen
/data/shard-msvcf/shard-88/subshard-1/postproc/pgen/dragen.pvar
/data/shard-msvcf/shard-88/subshard-1/postproc/pgen/dragen.psam
/data/shard-msvcf/shard-89/subshard-1/postproc/pgen/dragen.pgen
/data/shard-msvcf/shard-89/subshard-1/postproc/pgen/dragen.pvar
/data/shard-msvcf/shard-89/subshard-1/postproc/pgen/dragen.psam
```

### Output data of `step4`

The output data of `step4` is stored in the project `<region>-dragen-igg-<timestamp>-msvcf-version-<version>` under the folder `/data/chrom-msvcf`. The files are split by chromosome and format.

The chromosome-level full msVCFs and site VCFs are located in the following folders.

```bash
/data/chrom-msvcf/chrom-20/full-vcf/dragen.vcf.gz
/data/chrom-msvcf/chrom-20/full-vcf/dragen.sites.vcf.gz
```

The chromosome-level post-processed reduced msVCF files (biallelic with cohort-level ML filtering) and the site VCF files are located in the following folders.

```bash
/data/chrom-msvcf/chrom-20/postproc-vcf/dragen.vcf.gz
/data/chrom-msvcf/chrom-20/postproc-vcf/dragen.sites.vcf.gz
```

The chromosome-level post-processed PGEN files are located in the following folders.

```bash
/data/chrom-msvcf/chrom-20/postproc-pgen/dragen.pgen
/data/chrom-msvcf/chrom-20/postproc-pgen/dragen.pvar
/data/chrom-msvcf/chrom-20/postproc-pgen/dragen.psam
```

### Output data of `step5`

The output data of `step5` is stored in the project `<region>-dragen-igg-<timestamp>-msvcf-version-<version>` under the folder `/data/chrom-msvcf`. The files are split by chromosome and format.

The chromosome-level annotation files of biallelic site VCF are located in the following folders.

```bash
/data/chrom-msvcf/chrom-20/full-anno/dragen.anno.json.gz
/data/chrom-msvcf/chrom-20/full-anno/dragen.anno.vcf.gz
```

## Details about the IGG demo

### Batch

In this IGG demo, we have two batches of samples: `batch-1` with 4 samples, and `batch-2` with 3 samples. The gVCFs are shared via the PopGen CLI bundle in your region and linked to the release demo project (see Installation section above). The assignment of samples to batches is given in the following files. Each file contains the URI path of gVCF files for samples in a specific batch.

```bash
cat project-data/batch-gvcf/batch-1.gvcfs.csv
ica://<region>-popgen-cli-release-demo/data/giab/v4.3.13f/HG001-NA12878-40x.hard-filtered.gvcf.gz
ica://<region>-popgen-cli-release-demo/data/giab/v4.3.13f/HG002.novaseq.pcr-free.35x.hard-filtered.gvcf.gz
ica://<region>-popgen-cli-release-demo/data/giab/v4.3.13f/HG003.novaseq.pcr-free.35x.hard-filtered.gvcf.gz
ica://<region>-popgen-cli-release-demo/data/giab/v4.3.13f/HG004.novaseq.pcr-free.35x.hard-filtered.gvcf.gz

$ cat project-data/batch-gvcf/batch-2.gvcfs.csv
ica://<region>-popgen-cli-release-demo/data/giab/v4.3.13f/HG005.hard-filtered.gvcf.gz
ica://<region>-popgen-cli-release-demo/data/giab/v4.3.13f/HG006_39x.hard-filtered.gvcf.gz
ica://<region>-popgen-cli-release-demo/data/giab/v4.3.13f/HG007_39x.hard-filtered.gvcf.gz
```

*TIP: You may modify the batch files and provide the URI of gVCF files from your own samples. We recommend 1000 samples per batch, which is a good balance between too few samples per batch (overhead of jobs and cost) and too many samples per batch (poor parallelization and high CPU RAM requirement per job).*

### Version

We use `batch-1` to demonstrate the original batch(es) of samples that are available for aggregation, and `batch-2` to represent a new set of samples available for iterative aggregation (N+1 scenario). We create two versions (or snapshots) of aggregation. The first version (version-1), with only `batch-1`, and the second version (`version-2`) with `batch-1` and `batch-2`. The versions are defined in the following files. Each file contains the batch IDs in a specific version.

```bash
cat project-data/version-batch/version-1.batches.txt
1

cat project-data/version-batch/version-2.batches.txt
1
2
```

*TIP: When you have multiple batches, you can define the versions accordingly. In an N+1 scenario, we recommend iterative aggregation with a few thousand or more samples per version. Since `step2`, `step3`, `step4` and `step5` are per version, rerunning these steps for only a couple of samples can result in high overhead of compute and cost without significant updates to the variants at the cohort level.*

### Gender information

For samples in each version, the sample gender information is provided in the following files.

```bash
$ cat project-data/version-gender/version-1.genders.tsv
#FID	IID	SEX
HG001-NA12878-40x	HG001-NA12878-40x	2
HG002.novaseq.pcr-free.35x	HG002.novaseq.pcr-free.35x	1
HG003.novaseq.pcr-free.35x	HG003.novaseq.pcr-free.35x	1
HG004.novaseq.pcr-free.35x	HG004.novaseq.pcr-free.35x	2

$ cat project-data/version-gender/version-2.genders.tsv
#FID	IID	SEX
HG001-NA12878-40x	HG001-NA12878-40x	2
HG002.novaseq.pcr-free.35x	HG002.novaseq.pcr-free.35x	1
HG003.novaseq.pcr-free.35x	HG003.novaseq.pcr-free.35x	1
HG004.novaseq.pcr-free.35x	HG004.novaseq.pcr-free.35x	2
HG005	HG005	1
HG006_39x	HG006_39x	1
HG007_39x	HG007_39x	2
```
The gender file follows the Plink2 sample format, tab-delimited, with 3 columns:

- `FID`: Family ID (this can be the same as IID for unrelated samples)
- `IID`: Individual ID
- `SEX`: 1=male, 2=female, 0=unknown

This information is required for ML cohort level filtering and msVCF to PGEN conversion.

*TIP: Please replace gender files with the gender information of your own samples following the same format.*

### gVCF version

As DRAGEN variant calling accuracy improves over DRAGEN versions, we recommend aggregating gVCFs from the same DRAGEN versions (major and minor) to avoid batch effects in variant accuracy across different DRAGEN versions. This batch effect may introduce noise in downstream analysis like GWAS or phasing. Differences in DRAGEN patch versions, but with the same major and minor version, are acceptable.

When different versions of gVCFs are detected across different batches, the IGG analyses will fail.

The DRAGEN version of input gVCF is defined in the following file

```bash
$ cat project-data/config/gvcf-version.txt
4.3.13f
```

*TIP: Please update the gVCF version according to your own gVCF files. You may extract version information from the header line of one of the gVCF files, which contains `DRAGENCommandLine`. The version format is `<major>.<minor>.<patch>`. For instance, only `4.3.13f` from the following header line:*

```bash
##DRAGENCommandLine=<ID=dragen,Version="SW: 13.021.732.4.3.13f, HW: 13.021.732",Date="Fri Feb 07 19:22:30 UTC 2025",CommandLineOptions="--enable-map-align true --enable-map-align-output true --output-format CRAM --enable-duplicate-marking true --enable-variant-caller true --vc-emit-ref-confidence GVCF --vc-enable-vcf-output true --enable-cnv true --sample-sex female --cnv-segmentation-mode SLM --enable-sv true --repeat-genotype-enable true --repeat-genotype-use-catalog expanded --enable-vntr true --enable-hla true --enable-star-allele true --enable-targeted true --cnv-enable-self-normalization true --hla-enable-class-2 true --force --lic-instance-id-location /opt/instance-identity --intermediate-results-dir /scratch --logging-to-output-dir true --output-file-prefix HG001-NA12878-40x --output-directory HG001-NA12878-40x --fastq-list fastq-list.csv --fastq-list-sample-id HG001-NA12878-40x --ref-dir /scratch/hg38-alt_masked/DRAGEN/10">
```

### Shard

IGG analyses are parallelized by genomic regions, or shards. The region of each shard is defined in this file.

```bash
$ cat project-data/config/shards.txt
1	chr1:1-34599028
2	chr1:34599029-69198056
3	chr1:69198057-103797084
4	chr1:103797085-138396112
5	chr1:138396113-172995140
...
```

For the human genome hg38, we use 102 shards:

- shard-1 to shard-100: chr1-22,X,chrY
- shard-101: chrM
- shard-102: 3341 alt contigs

### Subshard

In order to maximize the parallelization of `step3` and, at the same time, dynamically scale with an increasingly large number of samples and variants, we use dynamic subsharding. Each `step3` analysis is done on a subshard of genomic region, and each subshard contains the same number of variants. As the number of samples to aggregate increases, the number of variants per subshard decreases, resulting in more parallelized analyses in `step3`. Compared to fixed-size subsharding, the dynamic subsharding strategy avoids variation in compute resource and runtime due to variable variant density in the genome.

For each shard, the number of subshards is determined during `step2`, based on the total number of variants called in the shard. When you run the following script after `step2` analyses are completed,

```bash
./3.submit.sh count_subshards $version_id $shard_ids
```

It downloads the subsharding information (the number of subshards per shard) from all `step2` analyses output and creates following shard/subshards file in the metadata folder, and upload to sync up with this metadata in ICA project.

```bash
$ cat project-data/version-subshard/version-1.subshards.csv
88,1
89,1
```

The definition of subshard regions can be found from the `step2` output folder in ICA project `<region>-dragen-igg-<timestamp>-msvcf-version-<version>`,

```bash
/data/global-census/shard-88/dragen.subshards.stats.tsv
/data/global-census/shard-89/dragen.subshards.stats.tsv
```

*TIP: As the number of subshards per shard changes dynamically, please use the script the sync-up the IGG metadata, and do not manually create or edit these files `project-data/version-subshard/version-<version>.subshards.csv`.*

### Chromosome

The mapping between the chromosome ID and chromosome name is defined in the following file

```
$ cat project-data/config/chroms.txt | head
1	chr1
2	chr2
3	chr3
4	chr4
...
21	chr21
22	chr22
23	chrX
24	chrY
25	chrM
26	chrZ
```

For the human genome, we define all ALT contigs that are not autosomes, sex chromosomes, or mitochondria as chrZ (chrom-26).

### Reference Genome

In the demo, the reference genome (hg38) is shared via the PopGen CLI demo which can be found in the release demo project `<region>-popgen-cli-release-demo`. The URI of the reference genome is specified in the following file

```bash
$ cat project-data/config/ref.csv
ica://<region>-popgen-cli-release-demo/data/ref/fasta/hg38.fa
```

*TIP: For non-human genomes, or other versions of the human reference genome, please provide the URI of the reference fasta in the above format, which you have access to in your ICA domain.*

### Annotation

In the demo, for the human genome, we provide Illumina-Connected-Annotation database files shared via the PopGen CLI demo which can be found in the release demo project `<region>-popgen-cli-release-demo`. The annotation database files and version are specified in the following file. Note that the first line with "#" is required, which defines the overall root folder of annotation database files in ICA.

```bash
$ cat project-data/config/annodb.txt
# ica://<region>-popgen-cli-release-demo/data/ica-nirvana-v3.26.0
ica://<region>-popgen-cli-release-demo/data/ica-nirvana-v3.26.0/Cache/GRCh38.Ensembl.ndb
ica://<region>-popgen-cli-release-demo/data/ica-nirvana-v3.26.0/Cache/GRCh38.RefSeq.ndb
ica://<region>-popgen-cli-release-demo/data/ica-nirvana-v3.26.0/Cache/GeneSymbols.ndb
ica://<region>-popgen-cli-release-demo/data/ica-nirvana-v3.26.0/References/Homo_sapiens.GRCh38.Nirvana.dat
...
```

*TIP: If you have your own customized annotation database files to add (which are compatible with Illumina-Connected-Annotation), you may add the URIs with the same format to the file `project/config/annodb.txt`.*

*TIP: For non-human genomes, if you do not want to run annotation, you may skip `step5` in the submission script `3.submit.sh`*

### Concatenation by output format

For demo purposes, we simplified the concatenation parallelization in `step4`. For large-scale aggregation, to fully parallelize the concatenation per output file type, you may modify the submission script `3.submit.sh` to submit one concatenation analysis per chromosome per output file type. Supported file types are:

- reduced-vcf
- reduced-site-vcf
- pgen
- vcf
- site-vcf

```bash
function step4 () {
    version_id=$1
    chrom_ids=$2
    step_name="concat-msvcf"
    for concat_options in reduced-site-vcf reduced-vcf pgen site-vcf vcf; do
        popgen-cli dragen-igg submit \
            --input-project-data-folder-path $input_project_data_folder_path \
            --input-project-config-file-path $input_project_config_file_path \
            --output-analysis-json-folder-path $output_analysis_json_folder_path \
            --step-name $step_name \
            --version-ids $version_id \
            --chrom-ids $chrom_ids \
            --concat-options $concat_options \
            --analysis-instance-scratch-space 524288
    done
}
```

*TIP: When you are aggregating large cohorts (>=100K samples), concatenation of full msVCF per chromosome will not scale, and the output will be hundreds of TB to a few PB per each msVCF. The htslib indexing will stop working (taking a few days or weeks). This is when the concept of one msVCF per chromosome stops making sense. It is therefore strongly recommended to use the subshard full msVCF as in the "/shard-msvcf" folder of project `<region>-dragen-igg-<timestamp>-msvcf-version-<version-id>`, without concatenation.*

*TIP: Adjust the analysis scratch space size according to the expected output size from concatenation analyses using the option `--analysis-instance-scratch-space`. The default value is 524288 (in MB). As a reference, the output of aggregation of reduced msVCF of 250K samples can be of size a few hundred GB (autosomes) to a few TB (chromosome X).*
