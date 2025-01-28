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
The next downstream three steps depend upon a samplesheet file with a list of sample names and flow cells. This is a headerless (no column names) TSV file with sample names in the first column and flow cells in the second. This can be prepared from the sequencing QC TSV output shown above like so:

Here is an example:
```
Chile_404	PAW33034
Chile_406	PAW73369
Chile_509	PAW61512
Chile_511	PAW71368
Chile_516	PAW71977
```
## Data organization per cohort
We inititially transfer raw ONT data to cohort folders (e.g., /data/CARDPB/data/RUSH) and then organize this data in the following manner (modified output of tree -d -L 2 in the cohort directory):
```
.
├── POD5
│   └── UAB_ADC006
│       └── PAY16825
├── SCRIPTS
└── SEQ_REPORTS
    └── UAB_ADC006
        └── PAY16825
```

Scripts for data processing and analysis are placed in the SCRIPTS subdirectory, while raw POD5s and sequencing reports from each run are moved into further sample name and flow cell ID nested subdirectories (e.g., POD5/UAB_ADC006/PAY16825). This is the directory structure used to guide downstream data processing steps (basecalling, mapping, and variant calling).

## Basecalling
We basecall raw ONT data in POD5 format to unmapped BAMs using the ONT basecaller dorado. We currently use version 0.9.0 with the R10.4.1 E8.2 400bps super-accurate basecalling model. We also call 5mC/5hmC modification in the process.

To ensure basecalling completes correctly, we compare the sizes of POD5s to corresponding unmapped BAMs to ensure they all fall along a line.
## Mapping to a human genome reference
We map all reads to the human genome reference GRCh38 using minimap2 before checking for sample swaps and calling variants.
## Checking for sample swaps (case specific)
Depending upon cohort, we check for sample swaps at both the initial flow cell and merged sample levels. The sample swap calling procedure depends upon first subsetting samples' mapped reads per chromosome, calling single nucleotide variants from these subsets in relatively fast manner with Clair3, concatenating variant calls per sample, merging sample variant calls with corresponding short read variant calls, and then creating a king table with Plink2.
# Variant calling
We then use the genome alignment mappings to call a number of different variants, including single nucleotide variants (SNVs), short tandem repeats (STRs), and structural variants (SVs), using methods further described below.
## Single nucleotide variants (SNVs)
## Short tandem repeats (STRs)
## Structural variants (SVs)
## Merging variant calls
## Annotating variant calls
## Compute resource use summary
| Data processing step | Dependencies | Memory | CPUs | GPUs | Local scratch | Time allocation | 
| -------------------- | ------------ | ------ | ---- | ---- | ------------- | --------------- |
| Basecalling | Dorado 0.9.0, pod5 0.3.6 | 120GB | 30 | 2 A100 | 50GB | 2 days |
| Mapping | Minimap2 2.28 | 120GB | 40 | n/a | | 1 day |
| SNV calling | DeepVariant 1.8.0 | | | | 
| STR calling | Vamos 2.1.5 | | | | 
| SV calling | Sniffles 2.5.3 | | | |
