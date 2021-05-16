# ngsr.annotations

Helper functions for processing genome annotations

## Installation

```r
# Install ngsr.annotations from GitHub
devtools::install_github("mkabza/ngsr.annotations")
```

## Usage

### Preparing reference genome data directory

  * `prepare_reference_ensembl()` prepares a reference genome data directory using Ensembl data
  * `prepare_reference_kallisto_rna_velocity()` adds files for kallisto RNA velocity analysis to the reference genome data directory

### Working with TxDb objects

  * `txdb_to_gtf()` exports a TxDb object to a clean GTF file
 

