PHAGE=config["Phage name"]
ID=config["ID"]
fasta=config["fasta file"]
# identifying, orienting and removing adapter sequences and poly-A tails using Pychopper and cutadapt.
rule pychopper:
    output: 
        expand("results/processed_fastq/pychopper/pychopper_{phage}_{ident}/{phage}_enriched_{ident}_full_length_output.fq", phage=PHAGE, ident=ID),
        expand("results/processed_fastq/pychopper/pychopper_{phage}_{ident}/{phage}_control_{ident}_full_length_output.fq", phage=PHAGE, ident=ID),
        expand("results/processed_fastq/pychopper/pychopper_{phage}_{ident}/report_enriched.pdf", phage=PHAGE, ident=ID),
        expand("results/processed_fastq/pychopper/pychopper_{phage}_{ident}/report_control.pdf", phage=PHAGE, ident=ID),
        expand("results/processed_fastq/pychopper/pychopper_{phage}_{ident}/statistics_enriched.tsv", phage=PHAGE, ident=ID),
        expand("results/processed_fastq/pychopper/pychopper_{phage}_{ident}/statistics_control.tsv", phage=PHAGE, ident=ID)
    params:
        enriched=config["enriched fastq"],
        control=config["control fastq"],
        phage=PHAGE,
        ident=ID
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        pychopper -r {output[2]} -S {output[4]} {params.enriched} {output[0]}

        pychopper -r {output[3]} -S {output[5]} {params.control} {output[1]}
        """

rule cutadapt:
    input: "results/processed_fastq/pychopper/pychopper_{phage}_{ident}/{phage}_{cond}_{ident}_full_length_output.fq"
    output: temp("results/processed_fastq/cutadapt/{phage}_{cond}_{ident}_cutadapt_2.fq")
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
            cutadapt -a A{{10}} -e 1 -j 0 -o {output} {input}
        """
rule cutadapt2:
    input: "results/processed_fastq/cutadapt/{phage}_{cond}_{ident}_cutadapt_2.fq"
    output: "results/processed_fastq/cutadapt/{phage}_{cond}_{ident}_cutadapt.fq"
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
            cutadapt -g TTTCTGTTGGTGCTGATATTGCTGGG -e 1 -j 0 -o {output} {input}
        """
# Mapping reads onto genome
rule minimap2:
    input: "results/processed_fastq/cutadapt/{phage}_{cond}_{ident}_cutadapt.fq"
    output: temp("results/alignments/{phage}_{cond}_{ident}_minimap.sam")
    conda:
        "../envs/env_read_mapping.yaml"
    params:
        fa=fasta
    shell:
        """
        minimap2 -ax map-ont -k14 -t 8 {params.fa} {input} > {output}
        """
# creates index file for fasta file when needed
rule fastaIndex:
    output: expand("{fa}.fai",fa=fasta)
    conda:
        "../envs/env_read_mapping.yaml"
    params:
        fa=fasta
    shell:
        """
        samtools faidx {params.fa}
        """
# Soft clipping and converting SAM file to sorted BAM file
rule clipping:
    input: "results/alignments/{phage}_{cond}_{ident}_minimap.sam", expand("{fa}.fai",fa=fasta)
    output: temp("results/alignments/{phage}_{cond}_{ident}_clipped.sam")
    conda:
        "../envs/env_read_mapping.yaml"
    params:
        fa=fasta
    shell:
        """
        samclip --max 10 --ref  {params.fa} < {input[0]} > {output}
        """
# convert SAM files to BAM files
rule samToBam:
    input: "results/alignments/{phage}_{cond}_{ident}_clipped.sam"
    output: temp("results/alignments/{phage}_{cond}_{ident}.bam")
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        samtools view -q 30 -S -b {input} > {output}
        """
# sort BAM files
rule sortBam:
    input: "results/alignments/{phage}_{cond}_{ident}.bam"
    output: "results/alignments/BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam"
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        samtools sort {input} -o {output}
        """
#index BAM files
rule indexSortedBAM:
    input: "results/alignments/BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam"
    output: "results/alignments/BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam.bai"
    conda:
        "../envs/env_read_mapping.yaml"
    shell:
        """
        samtools index {input}
        """