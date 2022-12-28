rule getGenomeFile:
    input: expand("{fa}.fai",fa=config["fasta file"])
    output: '{phage}_genome.txt'
    shell:
        """
            awk -v OFS='\t' {{'print $1,$2'}} {input} > {output}
        """
rule PosEffRatios:
    input: 
        '{phage}_peak_calling/{phage}_enriched_{ident}.3end.plus.peaks.oracle.narrowPeak.counts.clustered.csv',
        '{phage}_genome.txt',
        'BAM_files_{phage}/{phage}_enriched_{ident}.sorted.bam'
    output: 
        'TTS_{phage}_{ident}/eff_ratios_{phage}_{ident}.plus.drop.coverage'
    params:
        t=config["TTS threshold"]
    shell:
        """
        bedFile=$(bedtools bamtobed -i {input[2]} | sort -k1,1 -k2,2n)
        posBedFile=$(grep -w "+" <(echo "$bedFile"))
        [ ! -e {output} ] || rm {output}
        > {output}

        sort -t ',' -k15 -n -u {input[0]} | awk -F ',' 'NR>1 {{print $2, $15}}' | while read -r chr peak;
        do
            genomeCovFile=$(bedtools genomecov -g {input[1]} -i <(echo "$posBedFile") -d)

            upstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 < TTS && $2 >= TTS-20' <(echo "$genomeCovFile") | \
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS-20, total/NR}}')
            downstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 > TTS && $2 <= TTS+20' <(echo "$genomeCovFile") |\
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS+20, total/NR}}')
            upstreamValue=$(awk '{{ print $3 }}' <(echo "$upstreamTTSaverageCoverage"))

            if (( $(echo $upstreamValue'>'0 | bc -l) ));
            then
                paste -d " " <(echo "$upstreamTTSaverageCoverage") <(echo "$downstreamTTSaverageCoverage") |\
                awk -v chr=$chr '{{print chr, $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}}' |\
                awk '{{ if($9 >= {params.t}) {{ print }} }}' >> {output}
            fi
        done
        """
rule NegEffRatios:
    input: 
        '{phage}_peak_calling/{phage}_enriched_{ident}.3end.minus.peaks.oracle.narrowPeak.counts.clustered.csv',
        '{phage}_genome.txt',
        'BAM_files_{phage}/{phage}_enriched_{ident}.sorted.bam'
    output: 
        'TTS_{phage}_{ident}/eff_ratios_{phage}_{ident}.minus.drop.coverage'
    params:
        t=config["TTS threshold"]
    shell:
        """
        bedFile=$(bedtools bamtobed -i {input[2]} | sort -k1,1 -k2,2n)
        posBedFile=$(grep -w "-" <(echo "$bedFile"))
        [ ! -e {output} ] || rm {output}
        > {output}

        sort -t ',' -k15 -n -u {input[0]} | awk -F ',' 'NR>1 {{print $2, $15}}' | while read -r chr peak;
        do
            genomeCovFile=$(bedtools genomecov -g {input[1]} -i <(echo "$posBedFile") -d)

            downstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 < TTS && $2 >= TTS-20' <(echo "$genomeCovFile") | \
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS-20, total/NR}}')
            upstreamTTSaverageCoverage=$(awk -v TTS="$peak" '$2 > TTS && $2 <= TTS+20' <(echo "$genomeCovFile") |\
            awk -v TTS="$peak" '{{total += $3}} END {{print TTS, TTS+20, total/NR}}')
            upstreamValue=$(awk '{{ print $3 }}' <(echo "$upstreamTTSaverageCoverage"))

            if (( $(echo $upstreamValue'>'0 | bc -l) ));
            then
                paste -d " " <(echo "$upstreamTTSaverageCoverage") <(echo "$downstreamTTSaverageCoverage") |\
                awk -v chr=$chr '{{print chr, $1, $2, $3, $4, $5, $6, $6/$3, 1-$6/$3}}' |\
                awk '{{ if($9 >= {params.t}) {{ print }} }}' >> {output}
            fi
        done
        """
rule PosCoverageDropToBed:
    input: 'TTS_{phage}_{ident}/eff_ratios_{phage}_{ident}.plus.drop.coverage'
    output: 'TTS_{phage}_{ident}/TTS_{phage}_{ident}.plus.bed'
    params:
        up=config["TTS sequence extraction"]["upstream"],
        down=config["TTS sequence extraction"]["downstream"]
    shell:
        """
        awk -v up={params.up} -v down={params.down} -v FS='\t' -v OFS='\t' -F ' ' '{{print $1, $2 - up, $2 + down, "TTS_" NR, "+"}}' {input} | uniq > {output}
        sed -i 's/\"//g' {output}
        """
rule NegCoverageDropToBed:
    input: 'TTS_{phage}_{ident}/eff_ratios_{phage}_{ident}.minus.drop.coverage'
    output: 'TTS_{phage}_{ident}/TTS_{phage}_{ident}.minus.bed'
    params:
        up=config["TTS sequence extraction"]["upstream"],
        down=config["TTS sequence extraction"]["downstream"]
    shell:
        """
        awk -v up={params.up} -v down={params.down} -v FS='\t' -v OFS='\t' -F ' ' '{{print $1, $2 - down, $2 + up, "TTS_" NR, "-"}}' {input} | uniq > {output}
        sed -i 's/\"//g' {output}
        """
rule extractSequencesTTS:
    input: 
        'TTS_{phage}_{ident}/TTS_{phage}_{ident}.{sign}.bed'
    output:
        'TTS_{phage}_{ident}/TTS_seq_{phage}_{ident}.{sign}.fa.out'
    params:
        fasta=config["fasta file"]
    shell:
        """
        bedtools getfasta -fi {params.fasta} -bed {input} -fo {output} -s -name 
        """