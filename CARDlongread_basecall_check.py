#!/usr/bin/env python3
# script to check that ONT basecalling completed successfully
# compares sizes of input POD5s and output UBAMs and checks ratios with linear regression
# prints plot of regression line plus table of outliers
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from io import BytesIO
from io import StringIO
import sh
import argparse
import statsmodels.api as sm
# make regression line
# print scatter plot with regression line
# print table with outliers based on residuals
# set up command line argument parser
parser = argparse.ArgumentParser(description='This program checks for basecalling completeness by comparing the sizes of input POD5 to output UBAM files. It outputs a scatterplot of POD5 and UBAM sizes and a table of outlier runs based on linear regression residuals.')

# get input and output arguments
parser.add_argument('--pod5_path', action="store", type=str, help="Path to POD5 files/directories (to get POD5 per run directory sizes)")
parser.add_argument('--ubam_path', action="store", type=str, help="Path to UBAM files (to get UBAM per run directory sizes)")
parser.add_argument('--input_table', action="store", type=str, help="Input table with no header and following order of columns: POD5 size, POD5 name, UBAM size, UBAM name. Generate with du --block-size=1K.")
parser.add_argument('--output', action="store", type=str, default=None, dest="output", help="Filename prefix for output files (optional)")
parser.add_argument('--plot_title', action="store", type=str, default=None, dest="plot_title", help="Title for output plot (optional)")

# parse arguments
results = parser.parse_args()
# debugging output
# print(results.pod5_path)
# print(results.ubam_path)
# print(results.input_table)

# throw error if no input provided
if ((results.pod5_path is None) and (results.ubam_path is None) and (results.input_table is None)):
	quit("ERROR: No input path or table provided!")
# throw error if either pod5 path or ubam path not provided and input_table not provided
if (((results.pod5_path is None) or (results.ubam_path is None)) and (results.input_table is None)):
	quit("ERROR: Either POD5 (--pod5_path) or UBAM path (--ubam_path) wasn't provided!")
# warning message that paths are used instead of table if both paths and table provided
if (results.pod5_path and results.ubam_path and results.input_table) is not None:
    print("Both paths and input table provided. Using paths to get input/output sizes.")
    
# set default output filename
if results.output is None:
    results.output='output_basecall_check'

# control stucture determining how to make POD5/UBAM data frame for analysis
if (results.pod5_path and results.ubam_path) is not None:
    # read pod5 sizes into data frame
    pod5_sizes_string = sh.du("--block-size=1K", sh.glob(results.pod5_path + "/*/*"))
    pod5_sizes_df = pd.read_csv(StringIO(pod5_sizes_string),sep='\t',header=None)
    # read ubam sizes into data frame
    ubam_sizes_string = sh.du("--block-size=1K", sh.glob(results.ubam_path + "/*/*.bam"))
    ubam_sizes_df = pd.read_csv(StringIO(ubam_sizes_string),sep='\t',header=None)
    # compare pod5 and ubam size tables to make sure they have the same number of rows
    # if not identical, quit with the respective error
    if pod5_sizes_df.shape[0] != ubam_sizes_df.shape[0]:
        quit('Error: POD5 and UBAM size lists are different lengths! Try again.')
    # make combined data frame
    pod5_ubam_sizes_df = pd.concat((pod5_sizes_df,ubam_sizes_df),axis=1,ignore_index=True)
# what to do if input table but not pod5/ubam paths
elif results.input_table is not None:
    # read data frame from table variable
    pod5_ubam_sizes_df = pd.read_csv(results.input_table,sep="\t",header=None)
# rename columns to make their names meaningful
pod5_ubam_sizes_df = pod5_ubam_sizes_df.rename(columns={0: "POD5 size", 1: "POD5 name", 2: "UBAM size", 3: "UBAM name"})
# calculate POD5/UBAM ratio
pod5_ubam_sizes_df["POD5/UBAM ratio"] = pod5_ubam_sizes_df["POD5 size"]/pod5_ubam_sizes_df["UBAM size"]
# print scatterplot with regression line
fig, ax = plt.subplots()
ax = sb.regplot(data=pod5_ubam_sizes_df,x="POD5 size",y="UBAM size")
# set axis labels
ax.set(xlabel="POD5 size (kb)",ylabel="UBAM size (kb)")
if results.plot_title is not None:
    ax.set(title=results.plot_title)
fig.savefig(results.output + "_scatterplot.png", format='png', dpi=150, bbox_inches='tight')
# close figure
fig.clf()
# print residuals plot
fig, ax = plt.subplots()
ax = sb.residplot(data=pod5_ubam_sizes_df,x="POD5 size",y="UBAM size")
# set axis labels
ax.set(xlabel="POD5 size (kb)",ylabel="UBAM size residual (kb)")
if results.plot_title is not None:
    ax.set(title=results.plot_title)
fig.savefig(results.output + "_residualplot.png", format='png', dpi=150, bbox_inches='tight')
# close figure
fig.clf()
# get standardized residuals for filtering
# response variable
y = pod5_ubam_sizes_df['UBAM size']
# predictor variable
x = pod5_ubam_sizes_df['POD5 size']
# add constant to predictor variable
x = sm.add_constant(x)
# fit linear regression model
model = sm.OLS(y,x).fit()
# get influence instance
influence=model.get_influence()
# obtain standardized residuals
standardized_residuals = influence.resid_studentized_internal
pod5_ubam_sizes_df['Standardized Residual']=standardized_residuals
# print standardized residual plot
fig, ax = plt.subplots()
ax = sb.scatterplot(data=pod5_ubam_sizes_df,x="POD5 size",y="Standardized Residual")
# set axis labels
ax.set(xlabel="POD5 size (kb)",ylabel="UBAM standardized residual (z-score)")
if results.plot_title is not None:
    ax.set(title=results.plot_title)
# add horizontal cutoffs for z scores of -2 (red), 0 (gray), and 2 (red)
ax.axhline(y=-2,color="red")
ax.axhline(y=0,color="gray")
ax.axhline(y=2,color="red")
fig.savefig(results.output + "_stdresidualplot.png", format='png', dpi=150, bbox_inches='tight')
# close figure
fig.clf()
# print table with outliers based on residuals
# filter for std residual with absolute value above 2
pod5_ubam_sizes_df_outliers = pod5_ubam_sizes_df[abs(pod5_ubam_sizes_df['Standardized Residual']) > 2]
# save as TSV file
pod5_ubam_sizes_df_outliers.to_csv(results.output + "_outliers.tsv",sep='\t',index=False)
# script complete
quit()
