rule getRandomSequencesFromFastaForTSS:
    output:
        'results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}.fa.out'
    params:
        fasta=config["phage fasta file"],
        upperBound=config["TSS info"]["bp_downstream"],
        lowerBound=config["TSS info"]["bp_upstream"],
        nSeq=config["numberRandomSeq"]
    shell:
        """
        Nseq={params.nSeq}
        samtools faidx {params.fasta}
        genomeName=$(cut -f1 {params.fasta}.fai)
        Nbp=$(cut -f2 {params.fasta}.fai)
        upperBound=$((Nbp - {params.lowerBound}))
        randomNumbers=$(shuf -i {params.lowerBound}-$upperBound -n $Nseq)

        NseqHalfPos=$((Nseq/2))
        NseqHalfPos=$(echo $NseqHalfPos | bc)
        NseqHalfNeg=$((Nseq - NseqHalfPos))

        randomNumbersPos=$(head -n $NseqHalfPos <(echo "$randomNumbers") | sort -n)
        randomNumbersNeg=$(tail -n $NseqHalfNeg <(echo "$randomNumbers") | sort -n)

        bedPos=$(awk -F'\t' -v OFS='\t' -v name=$genomeName '{{print name, $1 - {params.lowerBound}, $1 + {params.upperBound}, "TSS_" NR "_pos", 0, "+"}}' <<<$randomNumbersPos)
        bedNeg=$(awk -F'\t' -v OFS='\t' -v name=$genomeName '{{print name, $1 - {params.lowerBound}, $1 + {params.upperBound}, "TSS_" NR "_neg", 0, "-"}}' <<<$randomNumbersNeg)

        combinedBedFile=$(cat <(echo "$bedPos") <(echo "$bedNeg"))
        bedtools getfasta -fi {params.fasta} -bed <(echo "$combinedBedFile") -fo {output} -s -name
        """
rule getRandomSequencesFromFastaForTTS:
    output:
        'results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}.fa.out'
    params:
        fasta=config["phage fasta file"],
        upperBound=config["TTS info"]["bp_downstream"],
        lowerBound=config["TTS info"]["bp_upstream"],
        nSeq=config["numberRandomSeq"]
    shell:
        """
        Nseq={params.nSeq}
        samtools faidx {params.fasta}
        genomeName=$(cut -f1 {params.fasta}.fai)
        Nbp=$(cut -f2 {params.fasta}.fai)
        upperBound=$((Nbp - {params.lowerBound}))
        randomNumbers=$(shuf -i {params.lowerBound}-$upperBound -n $Nseq)

        NseqHalfPos=$((Nseq/2))
        NseqHalfPos=$(echo $NseqHalfPos | bc)
        NseqHalfNeg=$((Nseq - NseqHalfPos))

        randomNumbersPos=$(head -n $NseqHalfPos <(echo "$randomNumbers") | sort -n)
        randomNumbersNeg=$(tail -n $NseqHalfNeg <(echo "$randomNumbers") | sort -n)

        bedPos=$(awk -F'\t' -v OFS='\t' -v name=$genomeName '{{print name, $1 - {params.lowerBound}, $1 + {params.upperBound}, "TTS_" NR "_pos", 0, "+"}}' <<<$randomNumbersPos)
        bedNeg=$(awk -F'\t' -v OFS='\t' -v name=$genomeName '{{print name, $1 - {params.lowerBound}, $1 + {params.upperBound}, "TTS_" NR "_neg", 0, "-"}}' <<<$randomNumbersNeg)

        combinedBedFile=$(cat <(echo "$bedPos") <(echo "$bedNeg"))
        bedtools getfasta -fi {params.fasta} -bed <(echo "$combinedBedFile") -fo {output} -s -name
        """


rule calculateStability:
    input: "results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}.fa.out", "results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}.fa.out"
    output: 
        'results/{phage}_{ident}/TSS_{phage}_{ident}_stability.csv', 'results/{phage}_{ident}/TTS_{phage}_{ident}_stability.csv', 
        'results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}_stability.csv', 'results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}_stability.csv'
    params:
        fastaTSS= config["TSS info"]["fasta"],
        fastaTTS=config["TTS info"]["fasta"],
        upstreamTSS=config["TSS info"]["bp_upstream"],
        downstreamTSS=config["TSS info"]["bp_downstream"],
        upstreamTTS=config["TTS info"]["bp_upstream"],
        downstreamTTS=config["TTS info"]["bp_downstream"],
        smoothingWindow=20
    script:
        '../scripts/calculateStability.py'


rule calculateBendability:
    input: "results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}.fa.out", "results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}.fa.out"
    output: 
        'results/{phage}_{ident}/TSS_{phage}_{ident}_bendability.csv', 'results/{phage}_{ident}/TTS_{phage}_{ident}_bendability.csv','results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}_bendability.csv', 'results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}_bendability.csv'
    params:
        fastaTSS= config["TSS info"]["fasta"],
        fastaTTS=config["TTS info"]["fasta"],
        upstreamTSS=config["TSS info"]["bp_upstream"],
        downstreamTSS=config["TSS info"]["bp_downstream"],
        upstreamTTS=config["TTS info"]["bp_upstream"],
        downstreamTTS=config["TTS info"]["bp_downstream"],
        smoothingWindow=20
    script:
        '../scripts/calculateBendability.py'

rule calculateCurvature:
    input: "results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}.fa.out", "results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}.fa.out"
    output: 
        'results/{phage}_{ident}/TSS_{phage}_{ident}_curvature.csv', 'results/{phage}_{ident}/TTS_{phage}_{ident}_curvature.csv', 'results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}_curvature.csv', 'results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}_curvature.csv'
    params:
        fastaTSS= config["TSS info"]["fasta"],
        fastaTTS=config["TTS info"]["fasta"],
        upstreamTSS=config["TSS info"]["bp_upstream"],
        downstreamTSS=config["TSS info"]["bp_downstream"],
        upstreamTTS=config["TTS info"]["bp_upstream"],
        downstreamTTS=config["TTS info"]["bp_downstream"],
        smoothingWindow=10
    script:
        '../scripts/calculateCurvature.py'
rule plotTSS:
    input: 'results/{phage}_{ident}/TSS_{phage}_{ident}_{structprop}.csv', 'results/{phage}_{ident}/TSS_random_seq_{phage}_{ident}_{structprop}.csv'
    output: 'results/{phage}_{ident}/TSS_{phage}_{ident}_{structprop}.png', 'results/{phage}_{ident}/avg_TSS_{phage}_{ident}_{structprop}.png'
    params:
        upper=config["TSS info"]["bp_upstream_profile"],
        lower=config["TSS info"]["bp_downstream_profile"]
    script:
        '../scripts/plotStructuralProperties.py'
rule plotTTS:
    input: 'results/{phage}_{ident}/TTS_{phage}_{ident}_{structprop}.csv', 'results/{phage}_{ident}/TTS_random_seq_{phage}_{ident}_{structprop}.csv'
    output: 'results/{phage}_{ident}/TTS_{phage}_{ident}_{structprop}.png', 'results/{phage}_{ident}/avg_TTS_{phage}_{ident}_{structprop}.png'
    params:
        upper=config["TTS info"]["bp_upstream_profile"],
        lower=config["TTS info"]["bp_downstream_profile"]
    script:
        '../scripts/plotStructuralProperties.py'

