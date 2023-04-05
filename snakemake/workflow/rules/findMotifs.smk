rule findTSSmotifsMEME:
    input: 'results/TSS_{phage}_{ident}/TSS_seq_{phage}_enriched_{ident}.fa.out'
    output: 'results/TSS_{phage}_{ident}/meme_out/meme.html' 
    params:
        nmotifs=config["nmotifsTSS"]
    shell:
        """
        meme -dna -mod zoops -nmotifs {params.nmotifs} -oc results/TSS_{wildcards.phage}_{wildcards.ident}/meme_out -time 6000 {input}
        """
rule findTTSmotifsMEME:
    input: 'results/TTS_{phage}_{ident}/TTS_seq_{phage}_{ident}.fa.out'
    output: 'results/TTS_{phage}_{ident}/meme_out/meme.html' 
    params:
        nmotifs=config["nmotifsTTS"]
    shell:
        """
        meme -dna -mod zoops -nmotifs {params.nmotifs} -oc results/TTS_{wildcards.phage}_{wildcards.ident}/meme_out -time 6000 {input}
        """