configfile:	config["path_to_config_file"]

rule all:
	input: p_val = config["output_folder"] + "/binokulars_output/results/p_values.txt"

#---------------------------------------------------------------------------------------------------
# STEP 1: Generate Pileup
rule pileup:
	input:
		bam=config['path_to_bam'],
		ref_genome=config["path_to_reference_genome"]
	params:
		path=config['output_folder'] + "/biscuit_output",
		outdir=config['output_folder'],
		chr_file=config['chromosomes_file']	
	output:
		config["output_folder"]+ "/biscuit_output/my_pileup.vcf.gz"
	shell:	
		"""
		cd {params.outdir}
		mkdir -p biscuit_output
		mkdir -p binpolish/assets
		mkdir -p chromhmm/input_files
		mkdir -p chromhmm/output_files
		mkdir -p binokulars_output
		mkdir -p tmp

		n_lines=$(wc -l < {params.chr_file})
		if [ $n_lines -eq 1 ]
		then
		        chr=$(awk '{{print $1}}' {params.chr_file})
		else
		        chr="none"
		fi
		if [ $chr != "none" ]
		then
			biscuit pileup -g $chr -o {params.path}/my_pileup.vcf {input.ref_genome} {input.bam}	
		else
			biscuit pileup -o {params.path}/my_pileup.vcf {input.ref_genome} {input.bam}
		fi	
		bgzip {params.path}/my_pileup.vcf
		tabix -p vcf {output}
		"""

#---------------------------------------------------------------------------------------------------
# STEP 2: Generate Methylation Information

rule meth_info:
	input:
		ref_genome=config["path_to_reference_genome"],
		vcf_path=config["output_folder"]+"/biscuit_output/my_pileup.vcf.gz"
	params:
		temp_folder=config["output_folder"]+"/tmp",
	output:
		config['output_folder']+ "/biscuit_output/meth_info.bed"
	shell:
		"""
		biscuit vcf2bed -t cg {input.vcf_path} > temp.bed
		biscuit mergecg {input.ref_genome} temp.bed | awk '{{print $1 "\t" $2 "\t" ($3-1) "\t" $4 "\t" $5}}' > {output}
		mv temp.bed {params.temp_folder}
		"""

#---------------------------------------------------------------------------------------------------
# STEP 3: Format Methylation File

rule format_meth_data:
	input:
		config["output_folder"]+"/biscuit_output/meth_info.bed"
	params:
		script_path=config["path_to_jrc_seeker"]+"/scripts/format_methylation.py"
	output:
		meth_data=config["output_folder"]+"/biscuit_output/methylation_data.bed",
		meth_track=config["output_folder"]+"/biscuit_output/methylated_track.bed",
		unmeth_track=config["output_folder"]+"/biscuit_output/unmethylated_track.bed"
	shell:
		"python {params.script_path} -i {input} -o {output.meth_data} -m {output.meth_track} -u {output.unmeth_track}"



CHROMS = [i.split()[0] for i in open(config['chromosomes_file']) if 'random' not in i and 'M' not in i and 'alt' not in i and 'K' not in i and 'G' not in i and 'g' not in i and 'hap' not in i]
STATES=config['n_states']
SAMPLE=config['sample_name']


#---------------------------------------------------------------------------------------------------
# STEP 4: Binarize Methylation Data (ChromHMM Preprocessing)

rule binarize:
	input:
		meth_data=config["output_folder"]+"/biscuit_output/methylation_data.bed"
	params:
		bin_sz=config['bin_size'],
		im_low=config['lower_im_methylation_bound'],
		im_up=config['upper_im_methylation_bound'],
		sample=config['sample_name'],
		das=config['data_assignment_setting'],
		temp_folder=config['output_folder']+"/tmp",
		k=config['k_threshold'],
		chrom_lengths=config['path_to_chrom_length_file'],
		chromhmm_input=config['output_folder'] + "/chromhmm/input_files",
		chromhmm_output=config['output_folder'] + "/chromhmm/output_files",
		chromhmm = config['chromhmm'],
		script_path=config["path_to_jrc_seeker"]+"/scripts/binarize.py"

	output:
		expand("{path}/{sample_name}_{chrom}_binary.txt", path=config["output_folder"] + "/chromhmm/input_files", sample_name=SAMPLE, chrom=CHROMS)
	shell:
		"""
		python {params.script_path} -i {input.meth_data} -b {params.bin_sz} -l {params.im_low} -u {params.im_up} -s {params.sample} -d {params.das} -t {params.temp_folder} -k {params.k} -c {params.chrom_lengths} -o {params.chromhmm_input}
		"""


#---------------------------------------------------------------------------------------------------
# STEP 5: ChromHMM LearnModel

rule learnmodel:
	input:
		expand("{path}/{sample}_{chrom}_binary.txt", path=config["output_folder"] + "/chromhmm/input_files",sample=SAMPLE, chrom=CHROMS)
	params:
		n_states=config['n_states'],
		sample=config['sample_name'],
		chrom_lengths=config['path_to_chrom_length_file'],
		chromhmm_input=config['output_folder'] + "/chromhmm/input_files",
		chromhmm_output=config['output_folder'] + "/chromhmm/output_files",
		chromhmm = config['chromhmm'],
		chromhmm_iterations = config['chromhmm_it']
			
	output:	
		segments=expand("{path}/{sample_name}_{n_states}_segments.bed", path=config["output_folder"] + "/chromhmm/output_files",sample_name=SAMPLE, n_states=STATES),
		emission=expand("{path}/emissions_{n_states}.txt", path=config["output_folder"] + "/chromhmm/output_files",n_states=STATES)
	shell:
		"""
		java -mx4000M -jar {params.chromhmm} LearnModel -noautoopen -s 1 -r {params.chromhmm_iterations} {params.chromhmm_input} {params.chromhmm_output} {params.n_states} {params.chrom_lengths}
		"""

#---------------------------------------------------------------------------------------------------
# STEP 6: Process Emission Matrix

rule label_states:
	input:
		expand("{path}/emissions_{n_states}.txt", path=config["output_folder"]+"/chromhmm/output_files", n_states=STATES)
	params:
		script_path=config["path_to_jrc_seeker"]+"/scripts/label_states.py"
	output:
		config["output_folder"]+"/binpolish/assets/state_labels.txt"
	shell:
		"""
		python {params.script_path} -i {input} -o {output}
		"""

#---------------------------------------------------------------------------------------------------
# STEP 7: BinPolish

rule binpolish:
	input:
		segments=expand("{path}/{sample_name}_{n_states}_segments.bed", path=config["output_folder"] + "/chromhmm/output_files",sample_name=SAMPLE, n_states=STATES),
		labels=config["output_folder"]+"/binpolish/assets/state_labels.txt"

	params:
		blacklist=config['path_to_blacklist'],
		mappability=config['path_to_mappability_file'],
		segment_min_sz=config["segment_min_sz"],
		binpolishflank=config['binpolishflank'],
		script_path=config["path_to_jrc_seeker"]+"/scripts/binpolish.R"
	output:
		bp=config["output_folder"]+"/binpolish/polished_segmentation.bed",
		reg=config["output_folder"]+"/binpolish/im_regions.txt"
	shell:
		"""
		Rscript {params.script_path} -s {input.segments} -l {input.labels} -b {params.blacklist} -m {params.mappability} -y {params.segment_min_sz} -i {params.binpolishflank} -o {output.bp} -r {output.reg}
		"""

#---------------------------------------------------------------------------------------------------
# STEP 8: Generate Single Fragment Epireads

rule epiread:
	input:
		pileup=config["output_folder"]+"/biscuit_output/my_pileup.vcf.gz"
	params:	
		path=config["output_folder"]+"/biscuit_output",
		sample=config["sample_name"],
		bam=config['path_to_bam'],
		ref_genome=config["path_to_reference_genome"],
		temp=config['output_folder']+"/tmp",
		chr_file=config['chromosomes_file']
	output:
		epiread=config["output_folder"]+"/biscuit_output/" + config["sample_name"] + ".epiread",
		snps=config["output_folder"]+"/biscuit_output/" + config["sample_name"] + "_snps.bed"
	shell:
		"""
		biscuit vcf2bed -t snp {input.pileup} > {output.snps}
		n_lines=$(wc -l < {params.chr_file})
		if [ $n_lines -eq 1 ]
		then
		        chr=$(awk '{{print $1}}' {params.chr_file})
		else
		        chr="none"
		fi
		if [ $chr != "none" ]
		then
			biscuit epiread -g $chr -o {params.path}/temp.epiread -O -B {output.snps} {params.ref_genome} {params.bam} > {params.path}/temp.epiread
		else
			biscuit epiread -o {params.path}/temp.epiread -O -B {output.snps} {params.ref_genome} {params.bam} > {params.path}/temp.epiread 
		fi	
		
		sort -T {params.temp} -k2,2 -k3,3n {params.path}/temp.epiread | awk 'BEGIN{{ qname="" ; rec="" }} 
			qname == $2 {{ print rec"\t"$5"\t"$6"\t"$7"\t"$8 ; qname="" }} 
			qname != $2 {{ qname=$2 ; rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 ; pair=$3}}' > {output.epiread}
		
		rm {params.path}/temp.epiread
		"""


#---------------------------------------------------------------------------------------------------
# STEP 9: Process Epireads


rule process_epiread:
	input:
		epiread=config["output_folder"]+"/biscuit_output/" + config["sample_name"] + ".epiread"
	params:
		sample_name=config['sample_name'],
		out_path=config["output_folder"] + "/binokulars_output"

	output:
		out=expand("{path}/CHR/{chrom}", path=config["output_folder"] + "/binokulars_output", chrom=CHROMS)
	shell:
		"""
		cd {params.out_path}
		awk '{{ print $1 "\t" $3 "\t" gsub("C","C",$4) "\t" gsub("T","T",$4) "\t" gsub("C","C",$8) "\t" gsub("T","T",$8)}}' {input.epiread} > CT_reads.txt
		sort -k1,1 -k2n CT_reads.txt > CT_reads_sorted.txt
		mkdir -p CHR
		cp CT_reads_sorted.txt ./CHR
		cd CHR
		awk -F '\t' '{{print>$1}}' CT_reads_sorted.txt
		rm CT_reads_sorted.txt
		cd ../
		rm CT_reads.txt
		rm CT_reads_sorted.txt
		"""


#---------------------------------------------------------------------------------------------------
# STEP 10: Binokulars


rule binokulars:
	input:
		regions=config["output_folder"]+"/binpolish/im_regions.txt",
		epiread=expand("{path}/CHR/{chrom}", path=config["output_folder"] + "/binokulars_output",chrom=CHROMS)
	params:
		chr=config["output_folder"] + "/binokulars_output/CHR",
		path=config["output_folder"] + "/binokulars_output",
		l=config["bin_size"],
		N=config["permutation_iterations"],
		R=config["seed_binokulars"],
		c=config["binokulars_cores"],
		f=config["flank_length"],
		info=config["path_to_config_file"],
		out=config["output_folder"],
		path_binokulars=config["path_to_jrc_seeker"] + "/scripts/binokulars"
	output:
		p_val = config["output_folder"] + "/binokulars_output/results/p_values.txt"
	shell:
		"""
		cd {params.path}
		rm -rf "results"
		{params.path_binokulars} -t {input.regions} -i {params.chr} -l {params.l} -N {params.N} -f {params.f} -R {params.R} -o "results" -c {params.c}	
		"""
