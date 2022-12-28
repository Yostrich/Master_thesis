# Snakemake pipeline
## dependencies (link + recommended way to install)
* pychopper 
  * https://github.com/epi2me-labs/pychopper \
  * `mamba install -c epi2melabs -c conda-forge -c bioconda "epi2melabs::pychopper"`
* cutadapt 
  * https://cutadapt.readthedocs.io/en/stable/index.html
  * `mamba install -c bioconda cutadapt`
* minimap2
  * https://github.com/lh3/minimap2#install
  * `mamba install -c bioconda minimap2`
* samtools
  * http://www.htslib.org/
  * `mamba install -c bioconda samtools`
* bedtools
  * https://bedtools.readthedocs.io/en/latest/
  * `mamba install -c bioconda bedtools`
* termseq-peaks
  * https://github.com/nichd-bspc/termseq-peaks
```
git clone https://github.com/NICHD-BSPC/termseq-peaks
cd termseq-peaks
mamba install -c bioconda --file requirements.txt
python setup.py install
```
* r-base
  * https://www.r-project.org/
  * `mamba install -c conda-forge r-base`
