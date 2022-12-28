PHAGE=config["Phage name"]
ID=config["ID"]
fasta=config["fasta file"]

rule pychopper:
    output: 
        expand("pychopper_{phage}_{ident}/{phage}_enriched_{ident}_full_length_output.fq", phage=PHAGE, ident=ID),
        expand("pychopper_{phage}_{ident}/{phage}_control_{ident}_full_length_output.fq", phage=PHAGE, ident=ID)
    params:
        enriched=config["enriched fastq"],
        control=config["control fastq"],
        phage=PHAGE,
        ident=ID
    shell:
        """
        pychopper -r pychopper_{params.phage}_{params.ident}/{params.phage}_enriched_{params.ident}_report.pdf {params.enriched} {output[0]}

        pychopper -r pychopper_{params.phage}_{params.ident}/{params.phage}_control_{params.ident}_report.pdf {params.control} {output[1]}
        """

rule cutadapt:
    input: "pychopper_{phage}_{ident}/{phage}_{cond}_{ident}_full_length_output.fq"
    output: "{phage}_{cond}_{ident}_cutadapt.fq"
    shell:
        """
            temp=$(cutadapt -a A{{10}} -e 1 -j 0 {input})
            cutadapt -g TTTCTGTTGGTGCTGATATTGCTGGG -e 1 -j 0 -o {output} <(echo "$temp")
        """
rule minimap2:
    input: "{phage}_{cond}_{ident}_cutadapt.fq"
    output: "{phage}_{cond}_{ident}_minimap.sam"
    params:
        fa=fasta
    shell:
        """
        minimap2 -ax map-ont -k14 -t 8 {params.fa} {input} > {output}
        """
rule fastaIndex:
    output: expand("{fa}.fai",fa=fasta)
    params:
        fa=fasta
    shell:
        """
        samtools faidx {params.fa}
        """
rule clipping:
    input: "{phage}_{cond}_{ident}_minimap.sam", expand("{fa}.fai",fa=fasta)
    output: temp("{phage}_{cond}_{ident}_clipped.sam")
    params:
        fa=fasta
    shell:
        """
        samclip --max 10 --ref  {params.fa} < {input[0]} > {output}
        """
rule samToBam:
    input: "{phage}_{cond}_{ident}_clipped.sam"
    output: temp("{phage}_{cond}_{ident}.bam")
    shell:
        """
        samtools view -q 30 -S -b {input} > {output}
        """
rule sortBam:
    input: "{phage}_{cond}_{ident}.bam"
    output: "BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam"
    shell:
        """
        samtools sort {input} -o {output}
        """
rule indexSortedBAM:
    input: "BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam"
    output: "BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam.bai"
    shell:
        """
        samtools index {input}
        """