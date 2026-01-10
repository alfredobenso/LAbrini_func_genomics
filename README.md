# LAbrini_func_genomics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the computational analysis pipeline for genome-wide functional genomics of *Clostridium autoethanogenum* during autotrophic growth on syngas versus heterotrophic growth on fructose. The analysis includes:

- **Genome-wide identification of transcriptionally regulated (DTR) genes and differentially translated (DTL) sites**
- **Genome-wide identification of transcription initiation sites (TISs) and transcription 3' ends (T3Es)**
- **Genome-wide analysis of protein-RNA interactions** contributing to autotrophy in *C. autoethanogenum*

## Citation

If you use this code or data, please cite:

```
[Authors et al.]. Functional genomics of Clostridium autoethanogenum during autotrophic growth.
Nature Communications (2026). DOI: [to be added upon publication]
```

## Installation

### Requirements

- **R** (version ≥ 4.0)
- **Python** (version ≥ 3.8)

### R Dependencies

Install required R packages:

```r
install.packages(c("tidyverse", "DESeq2", "edgeR", "GenomicRanges", 
                   "GenomicAlignments", "Rsamtools", "ggplot2"))
```

### Python Dependencies

Install required Python packages:

```bash
pip install -r requirements.txt
```

## Repository Structure

### ./scripts

Contains scripts to reproduce the results shown in the manuscript:

#### `operon_unif.r`

Analyses the read-coverage profile over the genomic region spanned by a gene (operon).

#### `tis.r`

Identifies transcription initiation sites by the filtering procedure consisting of:

- Site filtering by site coverage and distance from start (for TIS) or stop (for T3E) codon of CDS
- Site filtering by the consistency of 5'UTR read-coverage profile with TIS usage
- Site filtering by read-coverage profile over the genomic region spanned by gene (operon)
- Site clustering by physical proximity (described in the Methods section)

#### `t3e.r`

Identifies transcription termination sites by the filtering procedure consisting of:

- Site filtering by site coverage and distance from start (for TIS) or stop (for T3E) codon of CDS
- Site filtering by read-coverage profile over the genomic region spanned by gene (operon)
- Site clustering by physical proximity (described in the Methods section)

#### `dtr.definition.r`

Analyses the reproducibility and clustering of RNA-seq data across *C. autoethanogenum* syngas and fructose chemostats and carries out differential transcription analysis.

#### `dtl.definition.py`

Implements three complementary metrics to quantify translational regulation:

- TL-STATUS
- TL-PROFILE
- POLY-FRACTIONS

#### `rna-seq.vs.pol-seq.r`

Provides:

- Reproducibility analysis of polysome-seq data across *C. autoethanogenum* syngas and fructose chemostats
- Differentially translated (DTL) genes using a consensus approach (CONSENSUS) over the three metrics implemented in `dtl.definition.py` and differentially transcribed (DTR) genes
- Overlap between up- or down-regulated differentially transcribed genes (DTR) and differentially translated genes (DTL) between syngas and fructose chemostats

#### `catrapid.protein-RNA.int.r`

Displays the protein-RNA interactions of the RBPs predicted by PRONTO-TK with the 5' UTRs and 3' UTRs of DTRs and DTLs.

### ./data

Contains input data required by the scripts. See `./data/README.md` for download instructions for large sequencing files.

### ./results

Contains the results produced by the scripts and is organized in:

- **DTL**: Results on DTLs definition and analysis
- **DTR**: Results on DTRs definition and analysis
- **TIS-T3E**: Results on TISs and T3Es identification
- **protein-RNA**: Results on protein-RNA interactions contributing to autotrophy

## Usage

### Running the Analysis Pipeline

1. Download the required data files as described in `./data/README.md`

2. Run the scripts in the following order:

```bash
# Transcriptional analysis
Rscript scripts/operon_unif.r
Rscript scripts/tis.r
Rscript scripts/t3e.r
Rscript scripts/dtr.definition.r

# Translational analysis
python scripts/dtl.definition.py
Rscript scripts/rna-seq.vs.pol-seq.r

# Protein-RNA interactions
Rscript scripts/catrapid.protein-RNA.int.r
```

## External Tools

This analysis makes use of the following published tools:

- **PRONTO-TK**: RNA-binding protein prediction  
  DOI: [10.1093/nargab/lqae112](https://doi.org/10.1093/nargab/lqae112)
  
- **catRAPID omics v2**: Protein-RNA interaction prediction  
  DOI: [10.1093/nar/gkab393](https://doi.org/10.1093/nar/gkab393)

Relevant information to reproduce RBP prediction in *C. autoethanogenum* is available in the **RBP** folder.

## Data Availability

- **Code and processed results**: This GitHub repository
- **Raw sequencing data**: See `./data/README.md` for access instructions
- **Reference genome**: *C. autoethanogenum* DSM 10061 (NCBI Assembly: GCA_040166795.1)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please:
- Open an issue on this GitHub repository
- Contact: [your email here]

## Acknowledgments

We thank all contributors and funding agencies that supported this work.
