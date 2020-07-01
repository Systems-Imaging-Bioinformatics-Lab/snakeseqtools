
#Call peaks using MACS2
if config["control_present"] == "yes":

  rule all:
      input:
          expand(directory("{sample}_macs2_control_{control}"),zip,sample=config["treatments"],control=config["controls"])


  rule macs2_control:
      input:
        treatments = "treatments/{sample}.pruned.bam",
        controls = "controls/{control}.pruned.bam" 
      output:
         o2 = directory("{sample}_macs2_control_{control}")
      params:
         name = "{sample}"
      log:
    	   macs2 = "log/{sample}_control_{control}.macs2"
	   
      shell:
         "macs2 callpeak -t {input.treatments} -c {input.controls} -n {params.name} {config[macs_params]} --outdir {output.o2} {config[broad_peak]} &> {log.macs2};"


else: 
  rule all:
    input:
        expand("{sample}_macs2/{sample}.BLfiltered.narrowPeak", sample=config["treatments"])
        
  rule macs2_no_control:
    input:
      treatments = "treatments/{sample}.pruned.bam"
    output:
       narrowpeak = "{sample}_macs2/{sample}_peaks.narrowPeak"
    params:
       name = "{sample}",
       dir = "{sample}_macs2"
    log:
       macs2 = "log/{sample}.macs2"
    shell:
       "macs2 callpeak -t {input.treatments} -n {params.name} {config[macs_params]} --outdir {params.dir} {config[broad_peak]} &> {log.macs2};"
  rule blacklist_filter:
    input:
        narrowpeak = rules.macs2_no_control.output.narrowpeak
    output:
        narrowpeak = "{sample}_macs2/{sample}.BLfiltered.narrowPeak",
        summits = "{sample}_macs2/sample}.summits.BLfiltered.bed"
    params:
        narrowpeak = "{sample}_macs2/{sample}_peaks.narrowPeak",
        summits = "{sample}_macs2/{sample}_summits.bed"
    shell:
        "bedtools intersect -a {params.narrowpeak} -b {config[blacklist]} -v > {output.narrowpeak} ; "
        "bedtools intersect -a {params.summits} -b {config[blacklist]} -v > {output.summits} ; "



