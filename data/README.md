# Data Directory

This directory should contain the input data required to run the analysis scripts. Due to file size limitations on GitHub (>100 MB), the raw sequencing data files are hosted externally.

## Download Instructions

**All data files are available at the following Dropbox link:**

🔗 **https://www.dropbox.com/scl/fo/hjuvcbiida2lfrg1rezzk/AFRCXzf9Up_thLw7ztwpKdk?rlkey=vvapo42ntsz6xg2thpslxkltq&st=jqtkg6gq&dl=0**

## Expected Directory Structure

After downloading and extracting the data, your `./data` directory should have the following structure:

```
data/
├── Transcriptome_raw_data/
│   └── [BAM files for RNA-seq analysis]
├── Polysome_seq/
│   └── [Polysome-seq data files]
├── GCA_040166795.1/
│   └── [Reference genome and annotation files]
└── [Additional data files required by scripts]
```

## File Descriptions

### Transcriptome_raw_data/

Contains BAM files from RNA-seq experiments for *C. autoethanogenum* grown on syngas and fructose chemostats.

**Required for scripts:**
- `tis.r` - Transcription initiation site identification
- `t3e.r` - Transcription termination site identification
- `dtr.definition.r` - Differential transcription analysis

**File format:** BAM (Binary Alignment Map)  
**Total size:** ~[SIZE] GB

### Polysome_seq/

Contains polysome-seq data for translational profiling.

**Required for scripts:**
- `dtl.definition.py` - Differential translation analysis
- `rna-seq.vs.pol-seq.r` - Integrated transcription-translation analysis

### GCA_040166795.1/

Reference genome assembly and annotation for *C. autoethanogenum* DSM 10061.

**Contents:**
- Genome FASTA file
- Gene annotation (GFF/GTF)
- Additional reference files

**Source:** NCBI Assembly GCA_040166795.1

## Data Access

1. Click the Dropbox link above
2. Click "Download" button to download the entire folder as a ZIP archive
3. Extract the ZIP file contents into the `./data` directory of this repository
4. Verify that all subdirectories are present

**Note:** You do NOT need a Dropbox account to download the files. The link provides direct public access.

## Data Size Information

- **Total compressed size:** ~[SIZE] GB
- **Total uncompressed size:** ~[SIZE] GB
- **Recommended disk space:** [SIZE] GB (including space for results)

## Alternative Access

If you encounter issues accessing the Dropbox link, please:
- Open an issue on this GitHub repository
- Contact the corresponding author(s) directly

## Data Reuse and Citation

The data provided here is released under [SPECIFY DATA LICENSE, e.g., CC-BY-4.0].

When reusing this data, please cite:
```
[Authors et al.]. Functional genomics of Clostridium autoethanogenum during autotrophic growth.
Nature Communications (2026). DOI: [to be added]
```

## Technical Notes

- All BAM files are sorted and indexed
- Files are compatible with standard bioinformatics tools (samtools, HTSeq, etc.)
- Quality control metrics are included in the manuscript supplementary materials
