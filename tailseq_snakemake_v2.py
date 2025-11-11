include: "tailseq_config_v2.py"

rule all:
	input:
		expand("analysis/tail/{sample}_analyzer_dummy.txt", sample=SAMPLES)

rule del_1bp:
	input:
		r1 = "fastq/{sample}_L001_R1_001.fastq.gz",
		r2 = "fastq/{sample}_L001_R2_001.fastq.gz"
	output:
		R1 = "fastq/end/{sample}_R1_1bp.fastq.gz",
		R2 = "fastq/end/{sample}_R2_1bp.fastq.gz"
	shell:
		"cutadapt -u -1 -U -1 -o {output.R1} -p {output.R2} {input.r1} {input.r2}"

rule remove_adapter:
	input:
		r1 = "fastq/end/{sample}_R1_1bp.fastq.gz",
		r2 = "fastq/end/{sample}_R2_1bp.fastq.gz"
	output:
		R1 = "fastq/trim/{sample}_R1_trim.fastq.gz",
		R2 = "fastq/trim/{sample}_R2_trim.fastq.gz"
	shell:
		"cutadapt -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' \
		-A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' \
		-o {output.R1} -p {output.R2} {input.r1} {input.r2}"

rule run_demultiplexer:
	input:
		r1 = "fastq/trim/{sample}_R1_trim.fastq.gz",
		r2 = "fastq/trim/{sample}_R2_trim.fastq.gz"
	params:
		manifest = DEMULTIPLEX_MANIFEST,
		demultiplexer_code = DEMULTIPLEXER
	output:
		demultiplex_dummy="analysis/demultiplex/{sample}_demultiplex_dummy.txt"
	shell:
		"python {params.demultiplexer_code} \
		-r1 {input.r1} -r2 {input.r2} \
		-m {params.manifest} -s {wildcards.sample} \
		-o analysis/demultiplex/"

rule analyze_tail:
	input:
		demultiplex_dummy = "analysis/demultiplex/{sample}_demultiplex_dummy.txt"
	params:
		demultiplex_dir = "analysis/demultiplex/",
		snRNA_reference = RNA_REF,
		tail_analyzer_code = TAIL_ANALYZER
	output:
		analyzer_dummy="analysis/tail/{sample}_analyzer_dummy.txt"
	shell:
		"python {params.tail_analyzer_code} \
		-i {params.demultiplex_dir} -d {params.snRNA_reference} \
		-n {wildcards.sample} -o analysis/tail/"



