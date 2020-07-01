configfile: "/nfs/turbo/umms-ukarvind/tstephie/CBSR/rnaseq/snakeRNAseq/config/callpeak_config.yaml"


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
        expand("{sample}_macs2", sample=config["treatments"])
        
  rule macs2_no_control:
    input:
      treatments = "treatments/{sample}.pruned.bam"
    output:
       o2 = directory("{sample}_macs2")
    params:
       name = "{sample}"
    log:
       macs2 = "log/{sample}.macs2"
    shell:
       "macs2 callpeak -t {input.treatments} -n {params.name} {config[macs_params]} --outdir {output.o2} {config[broad_peak]} &> {log.macs2};"



