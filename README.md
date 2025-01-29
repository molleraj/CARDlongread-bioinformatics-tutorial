# NIA CARD long read sequencing group bioinformatics tutorial
Summarizing the foundational data management, processing, and analysis steps we do in the NIA CARD LRS group.
# Data management
## Data transfer from the PromethIONs
We move raw ONT runs from the PromethIONs to Biowulf using Globus. We have Globus endpoints set up on each PromethION using Globus Connect Personal (https://www.globus.org/globus-connect-personal).
## Data backup to Google Cloud storage
While basecalling raw POD5 ONT data on Biowulf (see later in the tutorial), we also archive this data to the coldest-storage archival bucket on Google Cloud Platform. We transfer this data to Google Cloud through the Google Cloud Platform endpoint NIH HPC provides on Globus.
## Raw data sequencing QC
We use a set of two Python scripts to collect critical quality control parameters for weekly and cohort-wide sequencing runs, such as read N50, per flow cell output, experiment output, and flow cell starting active pores. More details can be found here: https://github.com/molleraj/longread-report-parser

Below are sample commands for running the sequencing report parser on a set of JSONs to generate a summary table and the QC dashboard script to create a spreadsheet containing a number of tables and 1D/2D visualizations of key QC measures like read N50, run data output, and starting active pores (initial count of pores available for sequencing).
```
# Execute with file list of json reports (one per line):
python3 CARDlongread_extract_from_json.py --filelist example_json_reports.txt --output example_output.tsv

# Alternatively, execute on all json files within a directory
# (does not descend into subdirectories)
python3 CARDlongread_extract_from_json.py --json_dir /data/CARDPB/data/PPMI/SEQ_REPORTS/example_json_reports/ --output example_output.tsv

# Make sequencing QC analytics spreadsheet from above QC output table (example_output.tsv)
python3 CARDlongread_extract_summary_statistics.py -input example_output.tsv -output example_summary_spreadsheet.xlsx -platform_qc example_platform_qc.csv -plot_title "PPMI tutorial example"
```
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

To ensure basecalling completes correctly, we compare the sizes of POD5s to corresponding unmapped BAMs to ensure they all fall along a line. We have previously done this by making a list of POD5 input and UBAM output sizes in kilobytes and plotting the input against the output sizes in Excel. Currently we are developing a python script to both prepare input and output size lists and then to find outliers relative to the regression line.  
## Mapping to a human genome reference
We map all reads to the human genome reference GRCh38 using minimap2 before checking for sample swaps and calling variants.
To perform basic QC after mapping, we use the tool cramino from the package nanopack (https://github.com/wdecoster/nanopack), which is available as a module on Biowulf.
## Checking for sample swaps (case specific)
Depending upon cohort, we check for sample swaps at both the initial flow cell and merged sample levels through whole genome alignment and variant calling. The sample swap calling procedure depends upon first subsetting samples' mapped reads per chromosome, calling single nucleotide variants from these subsets in relatively fast manner with Clair3, concatenating variant calls per sample, merging sample variant calls with corresponding short read variant calls, and then creating a king table with Plink2.
# Variant calling
We then use the long read genome alignment mappings to call a number of different variants, including single nucleotide variants (SNVs), short tandem repeats (STRs), and structural variants (SVs), using methods further described below.
## Single nucleotide variants (SNVs)
We call single nucleotide variants from long read alignments with several different tools, including Clair3, DeepVariant, and PEPPER-Margin-DeepVariant (PMDV).
## Short tandem repeats (STRs)
We call short tandem repeats from long read alignments using Vamos.
## Structural variants (SVs)
We call structural variants from long read alignments using Sniffles.
## Merging variant calls
Variant call merging depends upon variant type.
## Annotating variant calls
Like merging, we annotate variant calls with different tools depending upon variant type.
## Compute resource use summary
| Data processing step | Dependencies | Memory | CPUs | GPUs | Local scratch | Time allocation | 
| -------------------- | ------------ | ------ | ---- | ---- | ------------- | --------------- |
| Basecalling | Dorado 0.9.0, pod5 0.3.6 | 120GB | 30 | 2 A100 | 50GB | 2 days |
| Mapping | Minimap2 2.28, samtools 1.21 | 120GB | 40 | n/a | n/a | 1 day |
| SNV calling | DeepVariant 1.8.0, singularity | 64GB | 16 | n/a | 50GB | 1 day |
| STR calling | Vamos 2.1.5, singularity | | | | 
| SV calling | Sniffles 2.5.3 | | | |
