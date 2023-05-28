rule createNotebook:
    input: 
        'TSS_{phage}_{ident}/TSS_seq_{phage}_enriched_{ident}.{sign}.fa.out',
        'TSS_{phage}_{ident}/random_seq_{phage}_{ident}.fa.out'

    output: 
        'TSS_{phage}_{ident}/structural_properties/{phage}_{ident}_curvature.{sign}.csv',
        'TSS_{phage}_{ident}/structural_properties/{phage}_{ident}_stability.{sign}.csv',
        'TSS_{phage}_{ident}/structural_properties/{phage}_{ident}_bendability.{sign}.csv',
        'TSS_{phage}_{ident}/structural_properties/random_seq_{phage}_{ident}_curvature.{sign}.csv',
        'TSS_{phage}_{ident}/structural_properties/random_seq_{phage}_{ident}_stability.{sign}.csv',
        'TSS_{phage}_{ident}/structural_properties/random_seq_{phage}_{ident}_bendability.{sign}.csv',

    log:
        notebook ="logs/notebooks/{phage}_{ident}.{sign}.ipynb"
    notebook: '../../notebooks/structuralProperties.py.ipynb'
rule getRandomSequencesFromFasta:
    output:
        '{TxS}_{phage}_{ident}/random_seq_{phage}_{ident}.fa.out'
    params:
        fasta=config["fasta file"],
        upperBound=lambda w: config["{} sequence extraction".format(w.TxS)]["downstream"],
        lowerBound=lambda w: config["{} sequence extraction".format(w.TxS)]["upstream"]
    shell:
        """
        Nseq=201
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

        bedPos=$(awk -F'\t' -v OFS='\t' -v name=$genomeName '{{print name, $1 - {params.lowerBound}, $1 + {params.upperBound}, "{wildcards.TxS}_" NR "_pos", "+"}}' <<<$randomNumbersPos)
        bedNeg=$(awk -F'\t' -v OFS='\t' -v name=$genomeName '{{print name, $1 - {params.lowerBound}, $1 + {params.upperBound}, "{wildcards.TxS}_" NR "_neg", "-"}}' <<<$randomNumbersNeg)

        combinedBedFile=$(cat <(echo "$bedPos") <(echo "$bedNeg"))
        bedtools getfasta -fi {params.fasta} -bed <(echo "$combinedBedFile") -fo {output} -s -name
        """
