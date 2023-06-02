signToSymbolPos={"higher": "+", "lower": "-"}
signToInverseSymbolPos={"higher": "-", "lower": "+"}

# Retrieves all peaks from the clustered peaks files for both enriched and control
rule getPeaks:
    input: 
        '{phage}_peak_calling/{phage}_enriched_{ident}.5end.{sign}.peaks.oracle.narrowPeak.counts.clustered.csv',
        '{phage}_peak_calling/{phage}_control_{ident}.5end.{sign}.peaks.oracle.narrowPeak.counts.clustered.csv'
    output: 
        temp('{phage}_{ident}_enrichedPeaks.{sign}.csv'), temp('{phage}_{ident}_controlPeaks.{sign}.csv')
    shell:
        """
        awk -F ',' '{{print $15}}' {input[0]}| awk 'NR>1' > {output[0]}
        awk -F ',' '{{print $15}}' {input[1]}| awk 'NR>1' > {output[1]}
        """
# Find all peaks that are found in both enriched and control and seperate them from non-overlapping peaks
rule splitOverlappingPeaks:
    input: '{phage}_{ident}_enrichedPeaks.{sign}.csv', '{phage}_{ident}_controlPeaks.{sign}.csv'
    output: temp('{phage}_{ident}_overlapping_peaks.{sign}.csv'), temp('{phage}_{ident}_remainingPeaks.{sign}.csv')
    shell:
        """
        grep -Fx -f {input[0]} {input[1]}|awk 'BEGIN{{FS="\t"; OFS=","}}$2=$1' > {output[0]}
        grep -Fxv -f {output[0]} {input[1]} > {output[1]}
        """
# Find peaks that do not exactly align in both conditions, but are within the error margin (PEAK-error:PEAK+error)
rule getOverlappingPeaksWithError:
    input: '{phage}_{ident}_remainingPeaks.{sign}.csv', '{phage}_{ident}_enrichedPeaks.{sign}.csv'
    output: temp('{phage}_{ident}_overlapping_peaks_{symb}_{error}.{sign}.csv')
    params:
        symb=lambda x: signToSymbolPos[x.symb],
        invSymb=lambda x: signToInverseSymbolPos[x.symb]
    shell:
        """
        errorPeaks=$(awk -F ',' '{{print $1{params.symb}{wildcards.error}}}' {input[0]})
        test=$(grep -Fx -f {input[1]} <(echo "$errorPeaks") || echo "")
        awk -F ',' 'BEGIN{{FS="\t"; OFS=","}}{{if ($1 > 0) {{print $1,$1{params.invSymb}{wildcards.error}}}}}' <(echo "$test") > {output}
        """
# Combine peaks that overlap exactly and peaks that do not overlap exactly, but within error margin
rule getAllOverlappingPeaks:
    input: '{phage}_{ident}_overlapping_peaks.{sign}.csv',
        expand('{{phage}}_{{ident}}_overlapping_peaks_{symb}_{error}.{{sign}}.csv', symb=["higher", "lower"], error=range(1,config["peak alignment error"]+1,1))
    output:
        temp('{phage}_{ident}_all_overlapping_peaks.{sign}.csv')
    shell:
        """
        cat {input}| awk 'NF' > {output}
        """
# Retrieves information about the peaks
rule getInformationFromOverlappingPeaks:
    input: 
        '{phage}_peak_calling/{phage}_enriched_{ident}.5end.{sign}.peaks.oracle.narrowPeak.counts.clustered.csv',
        '{phage}_peak_calling/{phage}_control_{ident}.5end.{sign}.peaks.oracle.narrowPeak.counts.clustered.csv',
        '{phage}_{ident}_all_overlapping_peaks.{sign}.csv'
    output:
        'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/enr_ratios_{phage}_{ident}.{sign}.csv'
    params:
        threshold=config["TSS Threshold"]
    shell:
        """
        commonEnriched=$(awk -F ',' '{{ print $1 }}' {input[2]} |\
        xargs -I {{}} awk -F "," '$15 == {{}} {{print $2, $3, $4, $15, $6, $18, $19}}' {input[0]})
        commonControl=$(awk -F ',' '{{ print $2 }}' {input[2]} |\
        xargs -I {{}} awk -F "," '$15 == {{}} {{print $2, $3, $4, $15, $6, $18, $19}}' {input[1]})

        paste -d ' ' <(echo "$commonEnriched") <(echo "$commonControl") |\
        awk '{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $7/$14}}' |\
        sort -g -k 2,2 | awk -v t={params.threshold} -F " " '$15 > t {{print }}' > {output}
        """
 # write information to BED files   
rule PosPeaksToBedFile:
    input: 'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/enr_ratios_{phage}_{ident}.plus.csv'
    output: 
        temp('results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_enriched_{ident}_peaks.plus.bed'),
        temp('results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_control_{ident}_peaks.plus.bed')
    params:
        up=config["TSS sequence extraction"]["upstream"],
        down=config["TSS sequence extraction"]["downstream"]
    shell:
        """
        awk -v FS='\t' -v OFS='\t' -F ' ' '{{print $1, $4 - {params.up}, $4 + {params.down}, "TSS_POS_" NR, 0, "+"}}' {input[0]} | uniq > {output[0]}
        awk -v FS='\t' -v OFS='\t' -F ' ' '{{print $8, $11 - {params.up}, $11 + {params.down}, "TSS_POS_" NR, 0, "+"}}' {input[0]} | uniq > {output[1]}
        sed -i 's/\"//g' {output[0]}
        sed -i 's/\"//g' {output[1]}
        """
rule NegPeaksToBedFile:
    input: 'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/enr_ratios_{phage}_{ident}.minus.csv'
    output: 
        temp('results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_enriched_{ident}_peaks.minus.bed'),
        temp('results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_control_{ident}_peaks.minus.bed')
    params:
        up=config["TSS sequence extraction"]["upstream"],
        down=config["TSS sequence extraction"]["downstream"]
    shell:
        """
        awk -v FS='\t' -v OFS='\t' -F ' ' '{{print $1, $4 - {params.down}, $4 + {params.up}, "TSS_NEG_" NR,0 , "-"}}' {input[0]} | uniq > {output[0]}
        awk -v FS='\t' -v OFS='\t' -F ' ' '{{print $8, $11 - {params.down}, $11 + {params.up}, "TSS_NEG_" NR,0 , "-"}}' {input[0]} | uniq > {output[1]}
        sed -i 's/\"//g' {output[0]}
        sed -i 's/\"//g' {output[1]}
        """
rule combinePosAndNegBedFiles:
    input: 'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_{cond}_{ident}_peaks.plus.bed', 'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_{cond}_{ident}_peaks.minus.bed'
    output: 'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_{cond}_{ident}_peaks.bed'
    shell:
        """
        cat {input} > {output}
        """
# extract sequences based on BED files
rule extractSequencesTSS:
    input: 
        'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/{phage}_{cond}_{ident}_peaks.bed'
    output:
        'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/TSS_seq_{phage}_{cond}_{ident}.fa.out'
    params:
        fasta=config["fasta file"]
    conda:
        "../envs/env_annotation.yaml"
    shell:
        """
        bedtools getfasta -fi {params.fasta} -bed {input} -fo {output} -s -name 
        """