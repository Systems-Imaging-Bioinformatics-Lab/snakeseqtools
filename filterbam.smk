configfile: "config/filter_config.yaml"


rule all:
  input:
    expand("pruned/{sample}.pruned.bam",sample=config["samples"])


if config["end_type"] == "pe" or config["end_type"] == "se":
  rule mark_duplicates:
    input:
        "bamfiles/{sample}.sorted.bam"
    output:
        bam = "mrkdup/{sample}.mrkdup.bam",
        metric = "mrkdup/{sample}.mrkdup.metric"
    shell:
        "export JAVA_OPTIONS=-Xmx12g ; "
        "PicardCommandLine MarkDuplicates I={input} O={output.bam} "
        "METRICS_FILE={output.metric} "
        "ASSUME_SORTED=True "
        "VALIDATION_STRINGENCY=LENIENT "

  rule index_dupmarked_bams:
    input:
        "mrkdup/{sample}.mrkdup.bam"
    output:
        "mrkdup/{sample}.mrkdup.bam.bai"
    shell:
        "samtools index {input} {output}"


  rule prune:
      input:
         bam = "mrkdup/{sample}.mrkdup.bam",
         bai = "mrkdup/{sample}.mrkdup.bam.bai"
      output:
         bam = "pruned/{sample}.pruned.bam",
         bai = "pruned/{sample}.pruned.bam.bai"
      params:
         flags = config['samtools_prune_flags']
      shell:
        "samtools view -b {params.flags} {input.bam} > {output.bam} ; "
        "samtools index {output.bam} {output.bai}"



