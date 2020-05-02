#SSvisualization.py - Scatter plots for overall summary set of unique substitutions
# Copyright (C) 2020  Evan Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib as mpl

from IPython.display import display

mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['savefig.dpi'] = 200

def read_overall_summary(config):
    overall_fpath = "{0}/summary/overall_summary.tsv".format(config['RUN']['RunName'])
    overall_df = pd.read_csv(overall_fpath,sep='\t',index_col=0)
    return overall_df

def get_scatter_alpha(overall_df):
    #Returns alpha value to use for scatter plots using overall_df
    alpha = 385 / len(overall_df)
    if alpha < 0.1:
        alpha *= 1.3
    return alpha

def AGS_13LGS_consensus(overall_df,outfile=""):
    """Returns DataFrame for which AGS Variant and Test Variant are equal. Assumes Test Variant is 13LGS.

    :param overall_df: DataFrame of run_wide unique substitutions
    :return: gs_consensus_df: DataFrame of subset of rows from overall_df where
    """
    gs_consensus_df = overall_df.loc[overall_df['Test Variant']==overall_df['AGS Variant'],:]
    if outfile:
        gs_consensus_df.to_csv(outfile,sep='\t',float_format='%.5f')
    return gs_consensus_df


def gene_split(gene_symbol, summary_df):
    #Splits summary_df by entries equal or not equal to gene_symbol
    gene_df = summary_df[summary_df["Gene"] == gene_symbol]
    rest_df = summary_df[~ (summary_df["Gene"] == gene_symbol)]
    return gene_df, rest_df

def scatter_setup(config,use_jsd_zscore=0,use_blos_zscore=0,title_str="",symbol="",
                  fname_identifier=""):
    """Generic parameters dictionary for setting up a scatter plot of unique substitutions JSD/ BLOSUM.

    :param config: see SSconfig
    :param use_jsd_zscore: 0,1,2 for raw score/ alignment Z-score/ unique substitution set-wide Z-score
    :param use_blos_zscore: 0,1,2 for raw score/ alignment Z-score/ unique substitution set-wide Z-score
    :param title_str: If not provided, results contains generic title string
    :param symbol: If provided, symbol is used to name figure file.
    :param fname_identifier: If provided, used to name figure file.
    :return: (dict) results: Contains information for title_str, x_field, y_field, figure_path
    """
    from SSutility.SSdirectory import create_directory
    run_name = config['RUN']['RunName']
    out_dir = "{0}/summary/plots".format(run_name)
    create_directory(out_dir)
    if not title_str:
        title_str = "JSD vs Average Test-Outgroup BLOSUM62 for all Unique Residues"
    plt.figure()

    # Mapping of parameters to column names for scatter plot data
    jsd_fields = {0: "JSD", 1: "JSD Alignment Z-Score", 2: "JSD US Z-Score"}
    blos_fields = {0: "Test-Outgroup BLOSUM62", 1: "Test-Outgroup BLOSUM Alignment Z-Score",
                   2: "Test-Outgroup BLOSUM US Z-Score"}
    x_field = jsd_fields[use_jsd_zscore]
    y_field = blos_fields[use_blos_zscore]
    # File Name adjustments
    jsd_fnames = {0: "JSD", 1: "JSD_alnZ", 2: "JSD_usZ"}
    blos_fnames = {0: "BLOS", 1: "BLOS_alnZ", 2: "BLOS_usZ"}
    path_suffix = "{0},{1}.png".format(jsd_fnames[use_jsd_zscore], blos_fnames[use_blos_zscore])
    if symbol:
        figure_path = "{0}/{1}_{2}".format(out_dir, symbol, path_suffix)
    elif fname_identifier:
        figure_path = "{0}/{1}_{2}".format(out_dir, fname_identifier, path_suffix)
    else:
        figure_path = "{0}/{1}".format(out_dir, path_suffix)

    labels = ['title_str','x_field','y_field','figure_path']
    params = dict(zip(labels,[title_str,x_field,y_field,figure_path]))
    return params

def generic_scatter(config,summary_df,use_jsd_zscore=0,use_blos_zscore=0,title_str="",dataset_id=""):
    """Generic scatter plot for summary_df.

    :param config:
    :param summary_df:
    :param (int) use_jsd_zscore: 0 for raw JSD, 1 for alignment specific Z-scores, 2 for unique dataset-wide Z-scores
    :param (int) use_blos_zscore: 0 for raw Test-Outgroup BLOSUM, 1 for alignment specific Z-scores,
    2 for unique dataset-wide Z-scores
    :param title_str: If provided, overwrites default figure title
    :param dataset_id: File name identifier for figure. Optional.
    :return: None, saves figure to path generated using figure_params
    """

    run_name = config['RUN']['RunName']
    figure_params = scatter_setup(config, use_jsd_zscore, use_blos_zscore, title_str,fname_identifier=dataset_id)
    labels = ['title_str', 'x_field', 'y_field', 'figure_path']
    title_str, x_field, y_field, figure_path = [figure_params[l] for l in labels]
    plt.scatter(x=summary_df[x_field], y=summary_df[y_field], alpha=alpha)
    plt.legend(["{0} Unique Substitutions".format(run_name)])
    plt.xlabel("Position " + x_field)
    plt.ylabel("Position {0}".format(y_field))
    plt.title(title_str)
    plt.savefig(figure_path)


def gene_split_scatter(config,symbol, summary_df, use_jsd_zscore=0,use_blos_zscore=0,
                       title_str=""):
    """Scatter plot with emphasized substitutions for summary entries corresponding to gene symbol.

    :param config: configparser object from config/config.txt
    :param symbol: gene symbol used to identify gene-specific data points for plot
    :param summary_df: DataFrame containing record
    :param (int) use_jsd_zscore: 0 for raw JSD, 1 for alignment specific Z-scores, 2 for unique dataset-wide Z-scores
    :param (int) use_blos_zscore: 0,1,2 similar to above for BLOSUM columns
    :param title_str: If provided, will be used as title string for plot.
    :return: None, saves figure to path generated using figure_params
    """
    run_name = config['RUN']['RunName']
    figure_params = scatter_setup(config,use_jsd_zscore,use_blos_zscore,title_str,symbol=symbol)
    labels = ['title_str', 'x_field', 'y_field', 'figure_path']
    title_str,x_field,y_field,figure_path = [figure_params[l] for l in labels]
    #Gene Split, plot each
    gene_df, rest_df = gene_split(symbol, summary_df)
    plt.scatter(x=rest_df[x_field], y=rest_df[y_field], alpha=alpha)
    plt.scatter(x=gene_df[x_field], y=gene_df[y_field], c="red", alpha=1)
    plt.legend([run_name,symbol])
    plt.xlabel("Position " + x_field)
    plt.ylabel("Position {0}".format(y_field))
    plt.title(title_str)
    plt.savefig(figure_path)


def summary_plots(config,plot_symbols=[]):
    # os.chdir("..")
    # from SSutility import config
    run_name = config['RUN']['RunName']
    gs_consensus_fpath = "{0}/summary/GS_consensus_uniques.tsv".format(run_name)

    overall_df = read_overall_summary(config)
    gs_consensus = AGS_13LGS_consensus(overall_df,outfile=gs_consensus_fpath)

    global alpha
    alpha = get_scatter_alpha(gs_consensus)
    global SCATTER_YLABEL
    SCATTER_YLABEL = "Average GS vs Outgroup BLOSUM62"

    generic_scatter(config,gs_consensus,use_jsd_zscore=0,use_blos_zscore=0,dataset_id="GS_consensus")
    generic_scatter(config, gs_consensus, use_jsd_zscore=2, use_blos_zscore=2, dataset_id="GS_consensus")


    for symbol in plot_symbols:
        gene_split_scatter(config,symbol,gs_consensus,use_jsd_zscore=0,use_blos_zscore=0)
        gene_split_scatter(config, symbol, gs_consensus, use_jsd_zscore=2, use_blos_zscore=2)
    # display(gs_consensus)

if __name__ == "__main__":
    # main()
    pass