configfile: "config.yaml"

rule all:
	input:
		expand("stringtie/mergedcounts/{sample}.gtf",sample=config["samples"]),
		expand("bws/SNR/{sample}.pos.SNR.bw", sample=config["samples"]),
        expand("bws/SNR/{sample}.neg.SNR.bw", sample=config["samples"]),
        expand("bws/paired/{sample}.pos.paired.bw", sample=config["samples"]),
        expand("bws/paired/{sample}.neg.paired.bw", sample=config["samples"])

rule bam_to_bigwig_subSno:
    input:
        bam="subSno/{sample}.subSno.bam",
        norm_factors="results/norm_factors.txt"
    output:
        pos_bw="bws/paired/{sample}.pos.paired.bw",
        neg_bw="bws/paired/{sample}.neg.paired.bw"
    threads: 8
    shell:
        """
        set -o pipefail
    	matched_line=$(awk -F'\t' -v sample="subSno.{wildcards.sample}.subSno.bam" '($1 == sample) {{print}}' {input.norm_factors})
    	echo "Matched line: $matched_line"
    	scaleFactor=$(awk -F'\t' -v sample="subSno.{wildcards.sample}.subSno.bam" '($1 == sample) {{print $2}}' {input.norm_factors})
    	echo "scaleFactor: $scaleFactor"
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.pos_bw} --filterRNAstrand forward --scaleFactor $scaleFactor -p {threads}
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.neg_bw} --filterRNAstrand reverse --scaleFactor $scaleFactor -p {threads}
        """
    
rule bam_to_bigwig_SNR:
    input:
        bam="SNR/{sample}.SNR_sorted.bam",
        norm_factors="results/norm_factors.txt"
    output:
        pos_bw="bws/SNR/{sample}.pos.SNR.bw",
        neg_bw="bws/SNR/{sample}.neg.SNR.bw"
    shell:
        """
        set -o pipefail
    	scaleFactor=$(awk -F'\t' -v sample="subSno.{wildcards.sample}.subSno.bam" '($1 == sample) {{print $2}}' {input.norm_factors})
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.pos_bw} --filterRNAstrand forward --scaleFactor $scaleFactor -p 8
        bamCoverage --bam {input.bam} -bs 1 --outFileFormat bigwig --outFileName {output.neg_bw} --filterRNAstrand reverse --scaleFactor $scaleFactor -p 8
        """

rule deseq2_normalization:
    input:
        counts="results/counts.featureCounts"
    output:
        "results/norm_factors.txt"
    shell:
        "Rscript scripts/deseq2_normalization.R {input.counts} {output}"

rule feature_counts:
    input:
        sam=expand("subSno/{sample}.subSno.bam", sample=config["samples"]), # list of sam or bam files
        annotation=config["annotation"],
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        multiext("results/{sample}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        1
    params:
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O -p --fracOverlap 0.2 -s 1"
    log:
        "logs/{sample}.log"
    wrapper:
        "0.72.0/bio/subread/featurecounts"

rule stringtie_merged_count:
	input:
		bam="subSno/{sample}.subSno.bam",
		gtf="stringtie/stringtiemerged.gtf"
	output:
		"stringtie/mergedcounts/{sample}.gtf"
	shell:
		"stringtie -e -B -p 8 -G {input.gtf} -o {output} {input.bam}"
        
rule stringtie_merge:
	input:
		gtf=expand("stringtie/{sample}.gtf",sample=config["samples"]),
		annotation=config["annotation"]
	output:
		"stringtie/stringtiemerged.gtf"
	shell:
		"stringtie --merge -i -p 8 -G {input.annotation} -o {output} {input.gtf}" 

rule stringtie:
	input:
		bam="subSno/{sample}.subSno.bam",
		annotation=config["annotation"]		
	output:
		"stringtie/{sample}.gtf"
	shell:
		"stringtie -p 8 -o {output} -G {input.annotation} --fr -m 100 {input.bam}"

rule getSNR:
	input:
		bam="subEE/{sample}.subEE.bam",
		bai="subEE/{sample}.subEE.bam.bai"
	output:
		"SNR/{sample}.SNR_sorted.bam",
		"SNR/{sample}.SNR_sorted.bam.bai"
	shell:
		"python scripts/get_SNR_bam.py -d ./SNR/ -s {wildcards.sample}.SNR -f {input.bam}"

rule indexsubEE:
    input:
        "subEE/{sample}.subEE.bam"
    output:
        temp("subEE/{sample}.subEE.bam.bai")
    log:
        "logs/samtools_index/{sample}.dedup.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"
		
rule subEE:
	input:
		bam="clipped/{sample}.clipped.bam",
		bai="clipped/{sample}.clipped.bam.bai"
	output:
		temp("subEE/{sample}.subEE.bam"),
	params:
		exonends_bed = config["exonends"]
	shell:
		"samtools view -L {params.exonends_bed} -U {output} -o /dev/null {input.bam}"
	
rule indexclipped:
    input:
        "clipped/{sample}.clipped.bam"
    output:
        temp("clipped/{sample}.clipped.bam.bai")
    log:
        "logs/samtools_index/{sample}.dedup.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"

rule removeclipping:
	input:
		bam="subSno/{sample}.subSno.bam",
		bai="subSno/{sample}.subSno.bam.bai"
	output:
		temp("clipped/{sample}.clipped.bam")
	shell:
		"python scripts/removeclipping.py {input.bam} {output}"
		
rule indexsubsno:
    input:
        "subSno/{sample}.subSno.bam"
    output:
        "subSno/{sample}.subSno.bam.bai"
    log:
        "logs/samtools_index/{sample}.dedup.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"

rule subsno:
	input:
		bam="dedup/{sample}.dedup.bam",
		bai="dedup/{sample}.dedup.bam.bai"
	output:
		"subSno/{sample}.subSno.bam"
	params:
		smallrnas_bed = config["smallRNAs"]
	shell:
		"samtools view -L {params.smallrnas_bed} -U {output} -o /dev/null {input.bam}"

rule samtools_index_dedup:
    input:
        "dedup/{sample}.dedup.bam"
    output:
        "dedup/{sample}.dedup.bam.bai"
    log:
        "logs/samtools_index/{sample}.dedup.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"

rule umitools_dedup:
	input:
		bam="aligned2/{sample}.sorted.bam",
		bai="aligned2/{sample}.sorted.bam.bai"
	output:
		"dedup/{sample}.dedup.bam"
	shell:
		"umi_tools dedup -I {input.bam} --paired --compresslevel 1 -S {output}"

rule samtools_index:
    input:
        "aligned2/{sample}.sorted.bam"
    output:
        "aligned2/{sample}.sorted.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@
    wrapper:
        "0.79.0/bio/samtools/index"

rule samtools_sort:
    input:
        "aligned2/{sample}.Unsort.bam"
    output:
        "aligned2/{sample}.sorted.bam"
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@.
    wrapper:
        "0.79.0/bio/samtools/sort"

rule samtools_view:
	input:
		 "aligned2/{sample}/Aligned.out.sam"
	output:
		temp("aligned2/{sample}.Unsort.bam")
	threads: 1
	shell:
		"samtools view -@ 8 -bhS {input} > {output}"

rule star_2pass:
    input:
        fq1 = "trimmed/{sample}.1.fastq.gz",
        fq2 = "trimmed/{sample}.2.fastq.gz",
        sjout = expand("aligned/{sample}/SJ.out.tab", sample=config["samples"]),
    output:
        temp("aligned2/{sample}/Aligned.out.sam")
    log:
        "logs/star/pe/{sample}.log"
    params:
        # path to STAR reference genome index
        index="/data/hyperion/ALBERT/A163B_Circadian/genome",
        # optional parameters
        extra="--alignMatesGapMax 100000 --outSAMstrandField intronMotif --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 1 --outSJfilterReads Unique -sjdbFileChrStartEnd {input.sjout}"
    threads: 4
    wrapper:
        "0.79.0/bio/star/align"

rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = "trimmed/{sample}.1.fastq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = "trimmed/{sample}.2.fastq.gz",
    output:
        temp("aligned/{sample}/Aligned.out.sam"),
        "aligned/{sample}/SJ.out.tab"
    log:
        "logs/star/pe/{sample}.log"
    params:
        # path to STAR reference genome index
        index=config["index"],
        # optional parameters
        extra="--alignMatesGapMax 100000 --outSAMstrandField intronMotif --outFilterMismatchNoverLmax 0.05 --outFilterMultimapNmax 1 --outSJfilterReads Unique"
    threads: 4
    wrapper:
        "0.79.0/bio/star/align"
    	
rule fastp_pe:
    input:
        sample=["umiextracted_reads/{sample}_R1.UMIextract.fq.gz", "umiextracted_reads/{sample}_R2.UMIextract.fq.gz"]
    output:
        trimmed=["trimmed/{sample}.1.fastq.gz", "trimmed/{sample}.2.fastq.gz"],
        # Unpaired reads separately
        unpaired1="trimmed/{sample}.u1.fastq.gz",
        unpaired2="trimmed/{sample}.u2.fastq.gz",
        html="report/pe/{sample}.html",
        json="report/pe/{sample}.json"
    log:
        "logs/fastp/pe/{sample}.log"
    params:
        adapters="--adapter_sequence AAGATCGGAAGAGC --adapter_sequence_r2 CTGTCTCTTATA",
        extra="--trim_poly_g --trim_poly_x -F 1"
    threads: 2
    wrapper:
        "v1.20.0/bio/fastp"

rule umiextract:
    input:
        read1=lambda wildcards: "{}_R1_001.fastq.gz".format(config["samples"][wildcards.sample], wildcards.sample),
        umi=lambda wildcards: "{}_R2_001.fastq.gz".format(config["samples"][wildcards.sample], wildcards.sample),
        read2=lambda wildcards: "{}_R3_001.fastq.gz".format(config["samples"][wildcards.sample], wildcards.sample),
    output:
        out1="umiextracted_reads/{sample}_R1.UMIextract.fq.gz",
        out2="umiextracted_reads/{sample}_R2.UMIextract.fq.gz",
    shell:
        "umi_tools extract -I {input.umi} --extract-method=regex --bc-pattern='^(?P<umi_1>.{{8}})' --read2-in={input.read1} --stdout=umiextracted_reads/filler.fq.gz --read2-out={output.out1} & umi_tools extract -I {input.umi} --extract-method=regex --bc-pattern='^(?P<umi_1>.{{8}})' --read2-in={input.read2} --stdout=/dev/null --read2-out={output.out2}"
