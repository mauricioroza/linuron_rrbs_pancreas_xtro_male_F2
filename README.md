# Pesticide-induced transgenerational alterations of genome-wide DNA methylation patterns in the pancreas of Xenopus tropicalis correlate with metabolic phenotypes

This repository contains the R code used in the analysis of RRBS data in the <a href="https://doi.org/10.1016/j.scitotenv.2024.170949" target="_blank">publication</a>

Raw RRBS fastq files have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number [PRJEB74774 ](https://www.ebi.ac.uk/ena/browser/view/PRJEB74774)

The sequences were processed using the [Nextflow nf-core/methylseq pipeline version 1.6.1](https://nf-co.re/methylseq/1.6.1)

The command used to launch the workflow was as follows:

```bash
nextflow run nf-core/methylseq -r 1.6.1 \
  --input '/linuron_f2_pancreas/*.fastq.gz' \
  --aligner bismark \
  --rrbs \
  --fasta GCF_000004195.4_UCB_Xtro_10.0_genomic.fasta \
  --single_end \
```

The Bismark .cov files were used as input, and the analysis was run in ./code/Lin_F2_Brain+Testis.rmd
