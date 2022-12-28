signToSymbol={"plus": "+", "minus": "-"}
numToTXS={"5": "TSS", "3": "TTS"}

rule createGenomecov:
    input: 'BAM_files_{phage}/{phage}_{cond}_{ident}.sorted.bam'
    output: '{phage}_peak_calling/{phage}_{cond}_{ident}.{num}end.{sign}.bedgraph'
    params:
        symb=lambda x: signToSymbol[x.sign]
    shell:
        """
        bedtools genomecov \
            -ibam {input} \
            -bga \
            -{wildcards.num} \
            -strand {params.symb} > {output[0]}
            """
rule termseqPeakCalling:
    input: '{phage}_peak_calling/{prefix}.{num}end.{sign}.bedgraph'
    output: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks'
    params:
        alpha=config["termseq alpha"],
        symb=lambda x: signToSymbol[x.sign]
    conda:
        '../envs/env_det_trans_bound.yaml'
    shell:
        """
        termseq_peaks {input} {input} --peaks {output} --strand {params.symb} -t {params.alpha}
        """
rule combineCovAndPeaks:
    input: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks', '{phage}_peak_calling/{prefix}.{num}end.{sign}.bedgraph'
    output: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts'
    shell:
        """
        bedtools intersect \
          -wao \
          -a {input[0]}.oracle.narrowPeak \
          -b {input[1]} \
          > {output}
        """
rule addTotalReadsToCounts:
    input: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts', 'BAM_files_{phage}/{prefix}.sorted.bam'
    output: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts.withTR'
    shell:
        """
        total_mapped=$(samtools view -c -F4 {input[1]})
        awk "{{print \$0, $total_mapped}}" {input[0]} > {output}
        """
rule clusterReads:
    input: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts.withTR'
    output: '{phage}_peak_calling/{prefix}.{num}end.{sign}.peaks.oracle.narrowPeak.counts.clustered.csv'
    params:
        cluster_width=lambda w: config["cluster width"]["{}".format(numToTXS[w.num])],
        minimum_coverage=lambda w: config["minimum coverage"]["{}".format(numToTXS[w.num])],
        symb=lambda x: signToSymbol[x.sign]
    shell:
        """
        chmod +x peak_clustering.r
        ./peak_clustering.r {input} {params.symb} {params.cluster_width} {params.minimum_coverage} {output}
        """

