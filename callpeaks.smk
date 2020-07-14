
#Call peaks using MACS2
if config["control_present"] == "yes":
	if config["broad_peak"] == "no":
	  rule all:
	      input:
	          expand("{sample}_macs2_control_{control}/{sample}.BLfiltered.narrowPeak",zip,sample=config["treatments"],control=config["controls"])

	  rule macs2_control:
	      input:
	        treatments = "treatments/{sample}.pruned.bam",
	        controls = "controls/{control}.pruned.bam" 
	      output:
	         narrowpeak = "{sample}_macs2_control_{control}/{sample}_peaks.narrowPeak"
	      params:
	         name = "{sample}",
	         dir = "{sample}_macs2_control_{control}"
	      log:
	    	   macs2 = "log/{sample}_control_{control}.macs2"
		   
	      shell:
	         "macs2 callpeak -t {input.treatments} -c {input.controls} -n {params.name} {config[macs_params]} --outdir {params.dir} &> {log.macs2};"

	  rule blacklist_filter:
		    input:
		        narrowpeak = rules.macs2_control.output.narrowpeak
		    output:
		        narrowpeak = "{sample}_macs2_control_{control}/{sample}.BLfiltered.narrowPeak",
		        summits = "{sample}_macs2_control_{control}/{sample}.summits.BLfiltered.bed"
		    params:
		        narrowpeak = "{sample}_macs2_control_{control}/{sample}_peaks.narrowPeak",
		        summits = "{sample}_macs2_control_{control}/{sample}_summits.bed"
		    shell:
		        "bedtools intersect -a {params.narrowpeak} -b {config[blacklist]} -v > {output.narrowpeak} ; "
		        "bedtools intersect -a {params.summits} -b {config[blacklist]} -v > {output.summits} ;"
	else:
	  rule all:
	      input:
	          expand("{sample}_macs2_control_{control}/{sample}.BLfiltered.broadPeak",zip,sample=config["treatments"],control=config["controls"])

	  rule macs2_control:
	      input:
	        treatments = "treatments/{sample}.pruned.bam",
	        controls = "controls/{control}.pruned.bam" 
	      output:
	         broadpeak = "{sample}_macs2_control_{control}/{sample}_peaks.broadPeak"
	      params:
	         name = "{sample}",
	         dir = "{sample}_macs2_control_{control}"
	      log:
	    	   macs2 = "log/{sample}_control_{control}.macs2"
		   
	      shell:
	         "macs2 callpeak -t {input.treatments} -c {input.controls} -n {params.name} {config[macs_params]} --outdir {params.dir} --broad &> {log.macs2};"

	  rule blacklist_filter:
		    input:
		        broadpeak = rules.macs2_control.output.broadpeak
		    output:
		        broadpeak = "{sample}_macs2_control_{control}/{sample}.BLfiltered.broadPeak"
		    params:
		        broadpeak = "{sample}_macs2_control_{control}/{sample}_peaks.broadPeak"
		    shell:
		        "bedtools intersect -a {params.broadpeak} -b {config[blacklist]} -v > {output.broadpeak}"		        


else:
	if config["broad_peak"] == "no":
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
		       "macs2 callpeak -t {input.treatments} -n {params.name} {config[macs_params]} --outdir {params.dir} &> {log.macs2};"
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
	else:
		  rule all:
		    input:
		        expand("{sample}_macs2/{sample}.BLfiltered.broadPeak", sample=config["treatments"])
		        
		  rule macs2_no_control:
		    input:
		      treatments = "treatments/{sample}.pruned.bam"
		    output:
		       broadpeak = "{sample}_macs2/{sample}_peaks.broadPeak"
		    params:
		       name = "{sample}",
		       dir = "{sample}_macs2"
		    log:
		       macs2 = "log/{sample}.macs2"
		    shell:
		       "macs2 callpeak -t {input.treatments} -n {params.name} {config[macs_params]} --outdir {params.dir} --broad &> {log.macs2};"
		  rule blacklist_filter:
		    input:
		        broadpeak = rules.macs2_no_control.output.broadpeak
		    output:
		        broadpeak = "{sample}_macs2/{sample}.BLfiltered.broadPeak"
		    params:
		        broadpeak = "{sample}_macs2/{sample}_peaks.broadPeak"
		    shell:
		        "bedtools intersect -a {params.broadpeak} -b {config[blacklist]} -v > {output.broadpeak} ; "
		       



