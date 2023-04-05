signToSymbol={"plus": "+", "minus": "-"}
numToTXS={"5": "TSS", "3": "TTS"}

# creates genome coverage file as input for the termseq peak calling
rule createGenomecov:
    input: 'results/BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam'
    output: temp('{phage}_peak_calling/{phage}_{cond}_{ident}.{num}end.{sign}.bedgraph')
    params:
        symb=lambda x: signToSymbol[x.sign]
    conda:
        "../envs/env_transcription_boundaries.yaml"
    shell:
        """
        bedtools genomecov \
            -ibam {input} \
            -bga \
            -{wildcards.num} \
            -strand {params.symb} > {output[0]}
            """
# peakcalling of the reads using termseqpeaks
rule termseqPeakCalling:
    input: '{phage}_peak_calling/{prefix}.{num}end.{sign}.bedgraph'
    output: temp('{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks'), temp('{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak')
    params:
        alpha=config["termseq alpha"],
        symb=lambda x: signToSymbol[x.sign]
    shell:
        """
        termseq_peaks {input} {input} --peaks {output[0]} --strand {params.symb} -t {params.alpha} --oracle-output {output[1]}
        """
# retrieves genome coverage information about the peaks found using termseqpeaks
rule combineCovAndPeaks:
    input: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak', '{phage}_peak_calling/{prefix}.{num}end.{sign}.bedgraph'
    output: temp('{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts')
    conda:
        "../envs/env_transcription_boundaries.yaml"
    shell:
        """
        bedtools intersect \
          -wao \
          -a {input[0]} \
          -b {input[1]} \
          > {output}
        """
# Adds total number of reads to the peaks coverage information.
rule addTotalReadsToCounts:
    input: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts', 'results/BAM_files_{phage}/{prefix}.sorted.bam'
    output: temp('{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts.withTR')
    conda:
        "../envs/env_transcription_boundaries.yaml"
    shell:
        """
        total_mapped=$(samtools view -c -F4 {input[1]})
        awk "{{print \$0, $total_mapped}}" {input[0]} > {output}
        """
# uses custom R script to cluster peaks based on peaks coverage information + total number of reads.
rule clusterReads:
    input: '{phage}_peak_calling/{phage}_{cond}_{ident}.{num}end.{sign}.peaks.oracle.narrowPeak.counts.withTR'
    output: '{phage}_peak_calling/{phage}_{cond}_{ident}.{num}end.{sign}.peaks.oracle.narrowPeak.counts.clustered.csv'
    conda:
        "../envs/env_transcription_boundaries.yaml"
    params:
        cluster_width=lambda w: config["cluster width"]["{}".format(numToTXS[w.num])],
        minimum_coverage=lambda w: config["minimum coverage"]["{}".format(w.cond)]["{}".format(numToTXS[w.num])],
        symb=lambda x: signToSymbol[x.sign]
    shell:
        """
        chmod +x peak_clustering.r
        ./peak_clustering.r {input} {params.symb} {params.cluster_width} {params.minimum_coverage} {output}
        """

