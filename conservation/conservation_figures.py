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
import seaborn as sns

from IPython.display import display
from utility.directoryUtility import taxid_dict,dir_vars

mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['savefig.dpi'] = 200

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

def generate_figures_dir():
    from utility import directoryUtility as dirutil
    summary_run_dir = dir_vars['summary_run']
    paths = ["{0}/figures/uniques_conservation"]
    for path in paths:
        subdir_path = path.format(summary_run_dir)
        for taxid in taxid_dict:
            tax_subpath = "{0}/{1}".format(subdir_path,taxid)
            dirutil.create_directory(tax_subpath)


def JSD_BLOS_scatter(taxid,jsd_a,blos_a,alpha,xlabel="JSD",ylabel="Normalized Test-Outgroup BLOSUM62",
                     jsd_threshs=[],blos_threshs=[],scatter_outpath=""):
    """Generic function for JSD/BLOSUM scatter plot. If threshs are provided, plot percentile lines. If scatter_outpath
    is provided, sace figure to path. """
    plt.figure()
    plt.xlim(0, 1)
    scatter = plt.scatter(x=jsd_a, y=blos_a, s=2, alpha=alpha)  # alpha)
    scatter.set_label("Non-gap uniques")
    if blos_threshs:
        xaxis_min, xaxis_max = 0, 1
        ax_80 = plt.hlines(blos_threshs[0], xaxis_min, xaxis_max, colors='r', linestyles='dotted', alpha=0.25)
        ax_95 = plt.hlines(blos_threshs[1], xaxis_min, xaxis_max, colors='r', linestyles='dashed', alpha=0.25)
        ax_80.set_label("80th percentile")
        ax_95.set_label("95th percentile")

    if jsd_threshs:
        yaxis_min, yaxis_max = plt.ylim()
        plt.ylim(yaxis_min, yaxis_max)
        ax_80 = plt.vlines(jsd_threshs[0], yaxis_min, yaxis_max, colors='r', linestyles='dotted', alpha=0.25)
        ax_95 = plt.vlines(jsd_threshs[1], yaxis_min, yaxis_max, colors='r', linestyles='dashed', alpha=0.25)
    plt.legend(loc="lower left", fontsize="x-small")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title("Non-gap unique substitutions: {0}".format(taxid_dict[taxid]))
    if scatter_outpath:
        plt.savefig(scatter_outpath)
    else:
        plt.show()

def metric_hist_plot(taxid, metric_label,metric_a, threshs=[],hist_fpath=""):
    plt.figure()
    plt.hist(metric_a,bins=40)
    plt.xlabel(metric_label)
    yaxis_min,yaxis_max = plt.ylim()
    ax_80 = plt.vlines(threshs[0], yaxis_min, yaxis_max, colors='r', linestyles='dotted', alpha=0.5)
    ax_95 = plt.vlines(threshs[1], yaxis_min, yaxis_max, colors='r', linestyles='dashed', alpha=0.5)
    ax_80.set_label("80th percentile")
    ax_95.set_label("95th percentile")
    plt.title("{0} Histogram: {1}".format(metric_label,taxid_dict[taxid]))
    plt.legend()
    if hist_fpath:
        plt.savefig(hist_fpath)
    else:
        plt.show()


def JSD_BLOS_dist_plots(taxid,overwrite=False):
    from conservation import conservation_summary as cs
    from conservation import analysis_calc as ac
    summary_run = dir_vars['summary_run']
    uc_dir = "{0}/figures/uniques_conservation/{1}".format(summary_run,taxid)
    print("Unique Conservation Dir: {0}".format(uc_dir))
    overall_tsv_fpath = "{0}/conservation/{1}_summary.tsv".format(summary_run,taxid)
    overall_df = cs.load_overall_summary_table(overall_tsv_fpath)

    non_gaps = cs.filter_gap_positions(overall_df)
    print("All uniques: {0}".format(len(overall_df)))
    print("Non_gaps: {0}".format(len(non_gaps)))
    ng_jsd= non_gaps['JSD']
    non_gaps_jsd_z_score = ac.calc_z_scores(non_gaps['JSD'])
    non_gaps_blos_z_score = ac.calc_z_scores(non_gaps['Test-Outgroup BLOSUM62'])
    ng_norm_blos = ((3-non_gaps['Test-Outgroup BLOSUM62'])/7)
    ng_norm_blos_nonna = ng_norm_blos.dropna()

    ng_jsd_blos = ng_jsd*ng_norm_blos
    ng_jsd_blos_nonna = ng_jsd_blos.dropna()
    ng_jsd_blos_threshs = [np.percentile(ng_jsd_blos_nonna,80),np.percentile(ng_jsd_blos_nonna,95)]
    alpha = 3850/(len(non_gaps))


    jsd_threshs = [np.percentile(ng_jsd,80),np.percentile(ng_jsd,95)]
    blos_threshs = [np.percentile(ng_norm_blos_nonna, 80), np.percentile(ng_norm_blos_nonna, 95)]
    print(jsd_threshs)
    print(blos_threshs)
    jsd_hist_fpath = "{0}/figures/uniques_conservation/{1}/{1}_jsd_hist.png".format(summary_run, taxid)
    if not os.path.exists(jsd_hist_fpath) or overwrite:
        metric_hist_plot(taxid,"JSD",ng_jsd,threshs=jsd_threshs,hist_fpath=jsd_hist_fpath)
    blos_hist_fpath = "{0}/figures/uniques_conservation/{1}/{1}_blos_hist.png".format(summary_run, taxid)
    if not os.path.exists(blos_hist_fpath) or overwrite:
        metric_hist_plot(taxid,"Normalized Test-Outgroup BLOSUM62",ng_norm_blos,threshs=blos_threshs,hist_fpath=blos_hist_fpath)
    jsd_blos_hist_fpath = "{0}/figures/uniques_conservation/{1}/{1}_jsd-blos_hist.png".format(summary_run, taxid)
    if not os.path.exists(jsd_blos_hist_fpath) or overwrite:
        metric_hist_plot(taxid, "Normalized JSD*BLOSUM62", ng_jsd_blos, threshs=ng_jsd_blos_threshs,
                         hist_fpath=jsd_blos_hist_fpath)

    # jsd_usz, blos_usz = overall_df['JSD US Z-Score'],overall_df['Test-Outgroup BLOSUM US Z-Score']
    # ax1 = sns.kdeplot(vals)
    scatter_outpath = "{0}/figures/uniques_conservation/{1}/{1}_norm_scatter.png".format(summary_run, taxid)
    if not os.path.exists(scatter_outpath) or overwrite:
        JSD_BLOS_scatter(taxid,ng_jsd,ng_norm_blos,alpha,jsd_threshs=jsd_threshs,blos_threshs=blos_threshs,
                         scatter_outpath=scatter_outpath)

def main():
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    generate_figures_dir()
    for taxid in taxid_dict:
        JSD_BLOS_dist_plots(taxid,overwrite=True)

if __name__ == "__main__":
    main()