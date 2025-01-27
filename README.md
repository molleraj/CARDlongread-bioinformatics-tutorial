# NIA CARD Long read sequencing group bioinformatics tutorial
Summarizing the foundational data management, processing, and analysis steps we do in the NIA CARD LRS group.
# Data management
## Data transfer from the PromethIONs
We move raw ONT runs from the PromethIONs to Biowulf using Globus. We have Globus endpoints set up on each 
## Data backup to Google Cloud storage
While basecalling
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
We basecall raw
## Mapping to a human genome reference
## Checking for sample swaps (case specific)
# Variant calling
## Single nucleotide variants (SNVs)
## Short tandem repeats (STRs)
## Structural variants (SVs)
## Merging variant calls
## Annotating variant calls
