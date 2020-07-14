# snakeseqtools
Analysis pipelines for RNAseq, ChIPseq and ATAC-seq


# Before running 
Modify all config files in the `config` folder to include sample details to be analyzed and required index paths as per specifications. You can also specify custom options for the trimming and alignment steps. Also, modify `sbat` files to set greatlakes account name and custom memory/core requirements.

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
