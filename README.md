# snakeseqtools
Analysis pipelines for RNAseq, ChIPseq and ATAC-seq. The steps are split up into 3 files which should be run in the following order if you start with fastq files:
1) `align.smk` : Runs QC, trim and align reads using desired method (bowtie2,hisat2 or kallisto)
2) `filterbam.smk`: Removes duplicates, and prunes bam files to remove unwanted reads using desired flags
3) `callpeaks.smk`: Calls peaks using MACS2 and removes blacklist regions

These scripts can be run as jobs in greatlakes using the corresponding `.sbat` files.

# Before running 
Modify all config files in the `config` folder to include sample details to be analyzed and required paths as per specifications. You can also specify custom options for the trimming and alignment steps. Also, modify `sbat` files to set greatlakes account name, email ID and custom memory/core requirements.

Samples can be run in single or paired ended modes, and the corresponding option can be specified in the `config.yaml` file as `se` or `pe` respectively. All input files must have the format `.fq.gz` and be placed a folder called `input` . Paired ended samples must be specified in the following format: `sample_1.fq.gz, sample_2.fq.gz`. Multiple samples can be specified under the keyword `samples` in the `config.yaml` file in the following format (note the space between `-` and `sample`):

```
samples:
  - sample1
  - sample2
  ...
```

Also, you can specify any one of the following alignment methods in the `config.yaml` file:
1) bowtie2
2) hisat2
3) kallisto

If you select `bowtie2` or `hisat2` the final output would be the sorted and indexed bam files. If you select `kallisto`, the final output would be the abundance files.

Before running `callpeaks.sbat`, create two folders named `treatments` and `controls` and place the pruned bam files from the previous steps in their respective folders. You can input multiple treatments and controls to call peaks by following the below mentioned `config` file format. You can also choose to run MACS2 without a control sample by updating the `callpeak_config.yaml` file appropriately. In this case, you do not have to create a `controls` folder.

If you have two treatment samples `t1.bam`, `t2.bam` and the corresponding controls are `c1.bam` and `c2.bam`, mention these in `callpeak_config.yaml` in the correct order:

```
treatments:
  - t1
  - t2
  ...
```

```
controls:
  - c1
  - c2
  ...
```

##### Notes on Filename Restrictions

Fastq files must be gzipped, and must have the extension `.fq.gz` for the pipelines to work.

To gzip all fastqs in a directory recursively, you can use the following command:

    Assuming fastq filenames end in ".fastq" here: If they end in ".fq", change the argument to -name "*.fq"
    
    find /path/to/fastq_dir -type f -name "*.fastq" -exec gzip {} \;

If your file extensions are not `.fq.gz`, for example if they are `.fastq.gz`, you can rename them all with the following command.

    find /path/to/fastq_dir -type f -name "*.fastq.gz" -exec rename .fastq.gz .fq.gz {} \;

### Dependiencies
* python 3
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [hisat2](http://daehwankimlab.github.io/hisat2/)
* [samtools/1.9](http://www.htslib.org/)



### How to run?

```bash
mkdir fastQC_output
Dry run: snakemake --snakefile <filename> --configfile config/<configfilename> -n
Actual run: sbatch <filename>.sbat
```
