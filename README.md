## indelicato_2023

### Skeletal muscle transcriptomics dissects the pathogenesis of Friedreich´s Ataxia (FRDA)

RNA sequencing data analysis of skeletal muscle biopsies from FRDA patients pre- and post-treatment with rhuEPO and from healthy controls.
Raw sequencing data can be downloaded from GEO ([GSE226646](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226646)) and should be placed in data/FRDA_01_RAW.
Alternatively, the the raw counts matrix can be downloaded for downstream analyses.


Conda environment for analyses:
```bash
conda env create -n frda-epo -f envs/frda-epo.yml
conda activate frda-epo
```

### Preprocessing

FastQC of the raw data
```bash
qsub preprocessing/fastqc.sh
```

Poly-A and quality trimming
```bash
source preprocessing/trimming.sh
```

Run the [nf-core/rnaseq](https://nf-co.re/rnaseq/3.6) pipeline (parameters and reference genome paths can be set in 'preprocessing/nf.config')
```bash
source preprocessing/preprocessing.sh
```

Run a check if paired samples are matched to the same donor with [NGSCheckMate-1.0.0](https://github.com/parklab/NGSCheckMate)
```bash
source preprocessing/sample_check.sh
```

### Analysis

Run the downstream analysis scripts
```bash
Rscript -e "rmarkdown::render('00_Preprocessing_Results.Rmd')"
Rscript -e "rmarkdown::render('01_Differential_Expression.Rmd')"
Rscript -e "rmarkdown::render('02_Pathway_Analysis.Rmd')"
Rscript -e "rmarkdown::render('03_Figures.Rmd')"
```

### Results

After execution of the analysis scripts, generated figures and tables can be found in 'results'.

