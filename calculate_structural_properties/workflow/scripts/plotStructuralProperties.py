import pandas as pd
import matplotlib.pyplot as plt

def plotStructuralProperty(inputRE, inputRandomSeq, output, lower, upper):
    structPropRE = pd.read_csv(inputRE, index_col=0)[upper:lower]
    structPropREavg = pd.read_csv(inputRE, index_col=0)[upper:lower].mean(axis=1)
    structPropRandomSeq = pd.read_csv(inputRandomSeq, index_col=0)[upper:lower].mean(axis=1)
    ax = structPropRE.plot(use_index=True, color="blue")
    ax = structPropRandomSeq.plot(use_index=True, color="black", ax=ax)
    ax.legend('')
    plt.savefig(output)
    plt.close()
    # ax2 = structPropRE.plot(use_index=True, color="blue")
    # ax2 = structPropRandomSeq.plot(use_index=True, color="black", ax=ax2)
    # ax2.legend('')
    # plt.savefig(outputavg)
def plotStructuralPropertyAverage(inputRE, inputRandomSeq, output, lower, upper):
    structPropRE = pd.read_csv(inputRE, index_col=0)[upper:lower].mean(axis=1)
    structPropRandomSeq = pd.read_csv(inputRandomSeq, index_col=0)[upper:lower].mean(axis=1)
    ax2 = structPropRE.plot(use_index=True, color="blue")
    ax2 = structPropRandomSeq.plot(use_index=True, color="black", ax=ax2)
    ax2.legend('')
    plt.savefig(output)

plotStructuralProperty(snakemake.input[0], snakemake.input[1], snakemake.output[0], snakemake.params["lower"], snakemake.params["upper"])
plotStructuralPropertyAverage(snakemake.input[0], snakemake.input[1], snakemake.output[1], snakemake.params["lower"], snakemake.params["upper"])