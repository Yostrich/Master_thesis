# Snakemake pipeline
## dependencies
The majority of dependencies are automatically installed in its own environments upon running the snakemake pipeline. Some dependencies, however need to be installed manually.
* snakemake
* termseq-peaks
  * https://github.com/nichd-bspc/termseq-peaks
```
git clone https://github.com/NICHD-BSPC/termseq-peaks
cd termseq-peaks
mamba install -c bioconda --file requirements.txt
python setup.py install
```
* MEME suite
  * https://meme-suite.org/meme/doc/install.html
