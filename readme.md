# Analysis code for Ailaoshan study

[![DOI](https://zenodo.org/badge/432829320.svg)](https://zenodo.org/badge/latestdoi/432829320)

## Citation

Yinqiu JI\*, Christopher CM BAKER\*,\*\*, Viorel D POPESCU, Jiaxin WANG, Chunying WU, Zhengyang WANG, Yuanheng LI, Lin WANG, Chaolang HUA, Zhongxing YANG, Chunyan YANG, Charles CY XU, Alex DIANA, Qingzhong WEN, Naomi E PIERCE\*\* and Douglas W YU\*\*. Measuring Protected-Area Effectiveness using Vertebrate Distributions from Leech iDNA. **Nature Communications**

\* - equal contributions<br/>
\*\* - corresponding authors

A preprint of this study is available from bioRxiv at [https://www.biorxiv.org/content/10.1101/2020.02.10.941336v3](https://www.biorxiv.org/content/10.1101/2020.02.10.941336v3) (doi:10.1101/2020.02.10.941336)

## Workflow and repository structure

The workflow for this study is documented using
[Snakemake](https://snakemake.github.io/). The [snakefile](snakefile) describes how results are generated. The Snakemake workflow is illustrated in the [DAG](docs/snakemake_dag.pdf), [filegraph](docs/snakemake_filegraph.pdf) and [rulegraph](docs/snakemake_rulegraph.pdf).

Note that some steps in this workflow are computationally intensive and/or have high memory requirements, especially model estimation and calculations on the posterior sample.

- [data](data) contains input data for this study as supplied by [Doug Yu](https://github.com/dougwyu). Note that [data/gis](data/gis) contains digital elevation model files downloaded from USGS Earth Explorer, and [data/pantheria](data/pantheria) contains data on mammals downloaded from the PanTHERIA database (see the [readme](data/pantheria/PanTHERIA_1-0_WR05_Aug2008_README.md) for details).

- [code](code) contains the R code for our analysis. [jags](jags) contains the JAGS code for the multi-species occupancy models, which we execute via the R package [jagsUI](https://cran.r-project.org/web/packages/jagsUI/index.html).

- [config](config) contains [config.yaml](config.yaml) which may be used to store user-specific parameters such as the API key used by [code/Ailaoshan_IUCNdata.R](code/Ailaoshan_IUCNdata.R) to access IUCN data.  An API key is not included in this repository, and you should supply your own if you want to use this code.

- Intermediate and final outputs get saved to [figures](figures), [tables](tables), [fasta](fasta), [rdata](rdata) and [rds](rds). In addition, [preOTU_networks](preOTU_networks) contains the correlation networks of preOTUs used by [Doug Yu](https://github.com/dougwyu) to generate the final OTUs. Some intermediate and final output files are included here to facilitate re-analysis of our results, but these can all be re-generated using the original input data and code in this repository.

- [scripts](scripts) and [slurm](slurm) may be used to locate files such as batch submission scripts, stdout and stderr files if running parts of the workflow in a cluster environment.

## Further information

- Please contact [Chris Baker](https://github.com/bakerccm) for any enquiries regarding this repository.
