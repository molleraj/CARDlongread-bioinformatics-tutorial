# NIA CARD Long read sequencing group bioinformatics tutorial
Summarizing the foundational data management, processing, and analysis steps we do in the NIA CARD LRS group.
# Data management
## Data transfer from the PromethIONs
We move raw ONT runs from the PromethIONs to Biowulf using Globus. We have Globus endpoints set up on each PromethION using Globus Connect Personal (https://www.globus.org/globus-connect-personal).
## Data backup to Google Cloud storage
While basecalling raw POD5 ONT data on Biowulf (see later in the tutorial), we also archive this data to the coldest-storage archival bucket on Google Cloud Platform. We transfer this data to Google Cloud through the Google Cloud Platform endpoint NIH HPC provides on Globus.
## Raw data sequencing QC
We use a set of two Python scripts to collect critical quality control parameters for weekly and cohort-wide sequencing runs, such as read N50, per flow cell output, experiment output, and flow cell starting active pores. More details can be found here: https://github.com/molleraj/longread-report-parser
# Data processing
## Sample sheet
The next downstream three steps depend upon a samplesheet file with a list of sample names and flow cells. This is a headerless (no column names) TSV file with sample names in the first column and flow cells in the second.

Here is an example:
```
Chile_404	PAW33034
Chile_406	PAW73369
Chile_509	PAW61512
Chile_511	PAW71368
Chile_516	PAW71977
```
## Data organization per cohort
## Basecalling
We basecall raw ONT data in POD5 format to unmapped BAMs using the ONT basecaller dorado. We currently use version 0.9.0 with the R10.4.1 E8.2 400bps super-accurate basecalling model. We also call 5mC/5hmC modification in the process.

To ensure basecalling completes correctly, we compare the sizes of POD5s to corresponding unmapped BAMs to ensure they all fall along a line.
## Mapping to a human genome reference
We map all reads to the human genome reference GRCh38 using minimap2 before checking for sample swaps and calling variants.
## Checking for sample swaps (case specific)
# Variant calling
## Single nucleotide variants (SNVs)
## Short tandem repeats (STRs)
## Structural variants (SVs)
## Merging variant calls
## Annotating variant calls
## Compute resource use summary
| Data processing step | Dependencies | Memory | CPUs | GPUs | Local scratch | 
| -------------------- | ------------ | ------ | ---- | ---- | ------------- |
| Basecalling | Dorado 0.9.0 | 120 GB | 30 | 2 A100 | 50GB |
