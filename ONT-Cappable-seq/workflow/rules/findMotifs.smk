rule findTSSmotifsMEME:
    input: 'results/transcript_boundaries/TSS_{phage}/TSS_{phage}_{ident}/TSS_seq_{phage}_enriched_{ident}.fa.out'
    output: 'results/MEME/TSS_{phage}_{ident}/meme_out/meme.html' 
    params:
        nmotifs=config["nmotifsTSS"]
    shell:
        """
        meme -dna -mod zoops -minw 25 -nmotifs {params.nmotifs} -oc results/TSS_{wildcards.phage}_{wildcards.ident}/meme_out -time 6000 {input}
        """
rule findTTSmotifsMEME:
    input: 'results/transcript_boundaries/TTS_{phage}/TTS_{phage}_{ident}/TTS_seq_{phage}_{ident}.fa.out'
    output: 'results/MEME/TTS_{phage}_{ident}/meme_out/meme.html' 
    params:
        nmotifs=config["nmotifsTTS"]
    shell:
        """
        meme -dna -mod zoops -minw 25 -nmotifs {params.nmotifs} -oc results/TTS_{wildcards.phage}_{wildcards.ident}/meme_out -time 6000 {input}
        """