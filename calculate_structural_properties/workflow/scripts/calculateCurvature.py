from dnacurve import CurvedDNA
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from operator import truediv
from scipy import signal

def getCurvature(seq):
    return CurvedDNA(seq, 'nucleosome').curvature[0,:].tolist()
def savePropertiesToCsv(propertyfunc, fastaFile, window,outCSV,bp_upstream,bp_downstream):
    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    df = pd.DataFrame(index=range(-1*bp_upstream,bp_downstream,1))
    for fasta in fasta_sequences:
        nameSeq, sequence = fasta.id, str(fasta.seq)
        df[str(nameSeq)] = propertyfunc(sequence)
        df[str(nameSeq)] = df[str(nameSeq)].rolling(window).mean()
    df.to_csv(outCSV, index=True)

savePropertiesToCsv(getCurvature, snakemake.params["fastaTSS"], snakemake.params["smoothingWindow"], snakemake.output[0], snakemake.params["upstreamTSS"], snakemake.params["downstreamTSS"])
savePropertiesToCsv(getCurvature, snakemake.params["fastaTTS"], snakemake.params["smoothingWindow"], snakemake.output[1], snakemake.params["upstreamTTS"], snakemake.params["downstreamTTS"])
savePropertiesToCsv(getCurvature, snakemake.input[0], snakemake.params["smoothingWindow"], snakemake.output[2], snakemake.params["upstreamTSS"], snakemake.params["downstreamTSS"])
savePropertiesToCsv(getCurvature, snakemake.input[1], snakemake.params["smoothingWindow"], snakemake.output[3], snakemake.params["upstreamTTS"], snakemake.params["downstreamTTS"])
