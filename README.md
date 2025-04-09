# NIA CARD long read sequencing group bioinformatics tutorial
Summarizing the foundational data management, processing, and analysis steps we do in the NIA CARD LRS group.
# Data management
## Data transfer from the PromethIONs
We move raw ONT runs from the PromethIONs to Biowulf using Globus. We have Globus endpoints set up on each PromethION using Globus Connect Personal (https://www.globus.org/globus-connect-personal).
## Data backup to Google Cloud storage
While basecalling raw POD5 ONT data on Biowulf (see later in the tutorial), we also archive this data to the coldest-storage archival bucket on Google Cloud Platform. We transfer this data to Google Cloud through the Google Cloud Platform endpoint NIH HPC provides on Globus.
## Raw data sequencing QC
We use a set of two Python scripts to collect critical quality control parameters for weekly and cohort-wide sequencing runs, such as read N50, per flow cell output, experiment output, and flow cell starting active pores. More details can be found here: https://github.com/molleraj/CARDlongread-report-parser

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
```
tail -n +2 example_output.tsv | cut -f1,6
```
The command above removes the output table header and then extracts column 1 and 6 of the output table (sample name and flow cell ID).

Here is an example sample sheet to be used with downstream analyses:
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

To ensure basecalling completes correctly, we compare the sizes of POD5s to corresponding unmapped BAMs to ensure they all fall along a line. We have previously done this by making a list of POD5 input and UBAM output sizes in kilobytes and plotting the input against the output sizes in Excel. More recently we developed a Python script to both prepare input and output size lists and then to find outliers relative to the regression line. This script is included in the repository (```CARDlongread_basecall_check.py```). It outputs a scatterplot of POD5 and UBAM sizes with a regression line, a residuals scatterplot, a standardized residuals scatterplot, and a table of outliers (|z-score| > 2 or user set cutoff) including POD5/UBAM ratio and standardized residual (z-score).

The Python script for checking basecalling completion through linear regression can be used like so:
```
usage: CARDlongread_basecall_check.py [-h] [--pod5_path POD5_PATH] [--ubam_path UBAM_PATH] [--input_table INPUT_TABLE] [--output OUTPUT] [--plot_title PLOT_TITLE] [--z_score_cutoff Z_SCORE_CUTOFF]

This program checks for basecalling completeness by comparing the sizes of input POD5 to output UBAM files. It outputs a scatterplot of POD5 and UBAM sizes and a table of outlier runs based on linear regression
residuals.

optional arguments:
  -h, --help            show this help message and exit
  --pod5_path POD5_PATH
                        Path to POD5 files/directories (to get POD5 per run directory sizes)
  --ubam_path UBAM_PATH
                        Path to UBAM files (to get UBAM per run directory sizes)
  --input_table INPUT_TABLE
                        Input table with no header and following order of columns: POD5 size, POD5 name, UBAM size, UBAM name. Generate with du --block-size=1K.
  --output OUTPUT       Filename prefix for output files (optional)
  --plot_title PLOT_TITLE
                        Title for output plot (optional)
  --z_score_cutoff Z_SCORE_CUTOFF
                        Z score cutoff for plots and reporting outliers (optional; default is 2)
```

Here are example scatterplots with the regression and standardized residuals:
<br></br>
<img src="https://github.com/user-attachments/assets/f9c6fc8e-4d6b-4fb0-8766-eb3eb3b99d91" width=720></img>
<br></br>
<img src="https://github.com/user-attachments/assets/a9e40ccc-c8af-4796-a100-61bebffac95b" width=720></img>
<br></br>
The red lines mark z-scores of 2 above and below the expected value, while the gray line marks a z-score of 0 (expected value based on the regression line).

## Mapping to a human genome reference
We map all reads to the human genome reference GRCh38 using minimap2 before checking for sample swaps and calling variants.
To perform basic QC after mapping, we use the tool cramino from the package nanopack (https://github.com/wdecoster/nanopack), which is available as a module on Biowulf. We have also developed a dashboard to aggregate cramino statistics over a group of alignments that we use for overall cohort QC (https://github.com/molleraj/CARDlongread-cramino-dashboard).
## Checking for sample swaps (case specific)
Depending upon cohort, we check for sample swaps at both the initial flow cell and merged sample levels through whole genome alignment and variant calling. The sample swap calling procedure depends upon first subsetting samples' mapped reads per chromosome, calling single nucleotide variants from these subsets in relatively fast manner with Clair3, concatenating variant calls per sample, merging sample variant calls with corresponding short read variant calls, and then creating a king table with Plink2.
# Variant calling
We then use the long read genome alignment mappings to call a number of different variants, including single nucleotide variants (SNVs), short tandem repeats (STRs), and structural variants (SVs), using methods further described below.
## Single nucleotide variants (SNVs)
We call single nucleotide variants from long read alignments with several different tools, including Clair3 (https://github.com/HKU-BAL/Clair3), DeepVariant (https://github.com/google/deepvariant), and PEPPER-Margin-DeepVariant (PMDV) (https://github.com/kishwarshafin/pepper). We have found through testing that either flow cell or merged sample (multiple flow cells combined) level alignment data should be subsetted by each chromosome in order to increase job throughput on the NIH HPC cluster and maintain reasonable resource (especially RAM) use on the nodes of the norm partition.
## Short tandem repeats (STRs)
We call short tandem repeats from long read alignments using Vamos (https://github.com/ChaissonLab/vamos). This is provided in a singularity container by the Chaisson lab.
## Structural variants (SVs)
We call structural variants from long read alignments using Sniffles (https://github.com/fritzsedlazeck/Sniffles). We first call variants in individual samples as Sniffles (SNF) output and then merge SNF files into a single VCF, also using Sniffles.
## Merging variant calls
Variant call merging depends upon variant type.
## Annotating variant calls
Like merging, we annotate variant calls with different tools depending upon variant type.
## Compute resource use summary
The estimates below are based on one 30x coverage sample (except for SNV calling - 1 sample alignment subsetted for a particular choromosome) and are generous especially concerning time (i.e., more resources provided than necessary).
| Data processing step | Dependencies | Memory | CPUs | GPUs | Local scratch | Time allocation | 
| -------------------- | ------------ | ------ | ---- | ---- | ------------- | --------------- |
| Basecalling | Dorado 0.9.0, pod5 0.3.6 | 120GB | 30 | 2 A100 | 50GB | 2 days |
| Mapping | Minimap2 2.28, samtools 1.21 | 120GB | 40 | n/a | n/a | 1 day |
| SNV calling | DeepVariant 1.8.0, singularity 4.1.5 | 64GB | 16 | n/a | 50GB | 1 day |
| STR calling | Vamos 2.1.5, singularity 4.1.5 | 32GB | 16 | n/a | n/a | 6 hours | 
| SV calling | Sniffles 2.5.3 | 16GB | 16 | n/a | n/a | 1 hour |
