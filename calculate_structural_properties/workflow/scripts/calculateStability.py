from dnacurve import CurvedDNA
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from operator import truediv
from scipy import signal

def seqTodiValues(seq, property_name):
    diNucTable = pd.read_csv('workflow/tables/dinucleotideTable.csv')
    value_list = []
    for i in range(0,len(seq)-1):
        j = i+2
        value_list.append(diNucTable.query(f"sequence == '{seq[i:j]}'")[property_name].item())
    return value_list

def savePropertiesToCsv(propertyfunc, fastaFile, window,outCSV,bp_upstream,bp_downstream):
    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    df = pd.DataFrame(index=range(-1*bp_upstream,bp_downstream-1,1))
    for fasta in fasta_sequences:
        nameSeq, sequence = fasta.id, str(fasta.seq)
        df[str(nameSeq)] = propertyfunc(sequence)
        df[str(nameSeq)] = df[str(nameSeq)].rolling(window).mean()
    df.to_csv(outCSV, index=True)

partial_stability = lambda y: seqTodiValues(y, "stability")
savePropertiesToCsv(partial_stability, snakemake.params["fastaTSS"], snakemake.params["smoothingWindow"], snakemake.output[0], snakemake.params["upstreamTSS"], snakemake.params["downstreamTSS"])
savePropertiesToCsv(partial_stability, snakemake.params["fastaTTS"], snakemake.params["smoothingWindow"], snakemake.output[1], snakemake.params["upstreamTTS"], snakemake.params["downstreamTTS"])
savePropertiesToCsv(partial_stability, snakemake.input[0], snakemake.params["smoothingWindow"], snakemake.output[2], snakemake.params["upstreamTSS"], snakemake.params["downstreamTSS"])
savePropertiesToCsv(partial_stability, snakemake.input[1], snakemake.params["smoothingWindow"], snakemake.output[3], snakemake.params["upstreamTTS"], snakemake.params["downstreamTTS"])
