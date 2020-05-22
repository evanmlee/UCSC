#analysis_calc.py - Metric calculations for variant occurence, BLOSUM scores,
#summary table writing and formatting
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
import numpy as np
import os

from IPython.display import display
from collections import OrderedDict
from conservation import jsd_calc

def gen_blos_df():
    from Bio.SubsMat.MatrixInfo import blosum62
    """Background distribution data from published code from below citation
    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    :return: aas: list of amino acid chars inlcuding gap char '-'
    :return: blosum62_bg: background distribution of amino acids from BLOSUM dataset
    :return: blos_df: DataFrame corresponding to BLOSUM62 matrix, index/cols are amino acid characters
    :return sim_matrix: np.ndarray of values in blos_df
    """
    aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','-']
    bg_probs = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, \
                0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
    blosum62_bg = bg_probs
    blos_df = pd.DataFrame(index=aas[:-1],columns=aas[:-1])
    for pair in blosum62:
        val = blosum62[pair]
        first, second = pair[0],pair[1]
        if first in aas and second in aas:
            blos_df.loc[first,second] = val
            blos_df.loc[second,first] = val
    sim_matrix = blos_df.values
    return aas, blosum62_bg, blos_df, sim_matrix

def indel_postiions(uniques,test_idx,unique_thresh):
    n_seq = len(uniques)
    indel_pos = []
    for pos in uniques:
        col = uniques.loc[:,pos]
        vc_metrics = variant_counts(col,test_idx)
        if vc_metrics['test_variant'] == '-' or vc_metrics['gap_fraction'] >= (1-unique_thresh/n_seq):
            indel_pos.append(pos)
    indel_cols = uniques.loc[:,indel_pos]
    return indel_cols

def id_exon_diffs(indel_cols, INDEL_THRESH=5):
    """Given a DataFrame of indel positions from align_df, identify exon_diffs as continuous stretches of at least
    INDEL_THRESH positions of adjacent indel positions.

    :param indel_cols: DataFrame of align_df columns corresponding to unique positions where test variant is '-' or
    outgroup variant is unviersally '-'
    :param INDEL_THRESH: int, threshold length for continuous stretches of indel positions to be returned.
    :return: exon_diffs, positions corresponding to continuous indel stretches as defined above.
    """
    # Iterate over all gap positions:
    # If current position already included in gap stretch positions, continue
    # Check if i+4 position corresponds to 5 sequential gaps starting at i
    # if so: continue checking sequential positions until non-sequential gap.
    # if not: move on to i+1 gap position
    #Side not i am so sorry if you are trying to read this function.
    indel_positions = indel_cols.columns
    exon_diffs = []
    for iter_idx in range(len(indel_positions) - 4):
        curr_gap = indel_positions[iter_idx]
        check_gap = indel_positions[iter_idx + 4]
        if curr_gap in exon_diffs:
            continue
        elif check_gap == curr_gap + 4:
            if iter_idx + 5 == len(indel_positions):
                prev_gap = check_gap
            else:
                for check_idx in range(iter_idx + 5, len(indel_positions)):
                    prev_gap = check_gap
                    check_gap = indel_positions[check_idx]
                    if check_gap == prev_gap + 1:
                        if check_idx == len(indel_positions) - 1:
                            prev_gap = check_gap
                        continue
                    else:
                        break
            for gap in range(curr_gap, prev_gap + 1):
                exon_diffs.append(gap)
        else:
            continue
    return exon_diffs

def filter_exon_diffs(uniques,test_species_idx,unique_thresh):
    """For a series of unique positions in align_df, split into exon differences (insertion or deletion stretches >=5
    aa) and unique substitutions (everything else).
    :param uniques:
    :param test_species_idx:
    :return:
    """
    indel_cols = indel_postiions(uniques,test_species_idx,unique_thresh)
    exon_diff_pos = id_exon_diffs(indel_cols)
    filt_uniques, exon_diffs = uniques.drop(columns=exon_diff_pos),uniques.loc[:,exon_diff_pos ]
    return filt_uniques, exon_diffs


def find_uniques(align_df, sub_freq_thresh, test_species_idx):
    """Identifies unnique positions from align_df with frequency <= sub_fre_threshold

    :param align_df: DataFrame of alignment characters. Columns: 1-indexed alignment positions, Index: record ids
    :param (int) sub_freq_thresh: max allowed number of instances of a substitution for it to be considered unique
    :param test_species_idx: record_id for test_species record in align_df for which uniques will be identified
    :param (boolean) display_uniques: If true, displays table of unique residues identified
    :return:
    """
    uniques = pd.DataFrame(index=align_df.index)
    for pos in align_df:
        col = align_df[pos]
        sub_counts = col.value_counts()
        test_var = col[test_species_idx][0]
        test_vc = sub_counts[test_var]
        if test_vc <= sub_freq_thresh and test_var in UNAMBIG_CHARS:
            uniques.loc[:,pos] = col
    return uniques

def gap_fraction(col):
    #Fraction of positions in col corresponding to gap character '-'
    return sum([c == '-' for c in col]) / len(col)

def sub_fraction(col, sub):
    #Fraction of positions in col corresponding to character sub.
    return sum([c == sub for c in col]) / len(col)

def align_pos_to_native_pos(align_df, idx, align_pos):
    """Converts alignment position to native sequence position (ie ignores gaps).

    :param align_df: DataFrame of msa characters
    :param idx: Index (row) of align_df to use for native position
    :param align_pos: alignment column position to convert to native position
    :return: native_pos: converted alignment position
    """
    pre_pos = align_df.loc[idx.values[0],:align_pos-1]
    pre_pos_vc = pre_pos.value_counts()
    native_pos = align_pos
    if '-' in pre_pos_vc:
        native_pos -= pre_pos_vc['-']
    return native_pos


def calc_z_scores(scores):
    """For an array/pandas Series of scores, calculates z-scores and returns them. Ignores NaN values.

    :param scores: array/Series of scores for which z-scores will be calculated
    :return: z_scores: array/Series of same dimensions/index as scores.
    """
    mean = np.nanmean(scores)
    std = np.nanstd(scores)
    z_scores = (scores-mean)/std
    return z_scores


def generate_jsd_series(test_idx,align_df,
                        keep_test_spec=False,use_gap_penalty=True):
    """ Generate a Jensen Shannon Difference Series for position indices within the MSA. JSD is by default based on
    align_df with test_species dropped (ie outgroup species). See JSDcalc.py for calculation details. Calculations
    based on below citation.

    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    https://compbio.cs.princeton.edu/conservation/
    :param test_species: index value in align_df corresponding to test_species. Row will be removed from calculation
    unless keep_test_spec set to True.
    :param align_df: DataFrame containing multiple sequence alignment (1-indexed positions as columns, record ids as index
    :param (boolean) keep_test_spec: Determines whether JSD calculation is done on outgroup or all records in align_df
    :return jsd_srs: pandas Series containing JSD values at each position (1-indexed)
    :return jsd_zscores: pandas Series containing JSD z-score (calculated using mean and std of this alignment only)
    """
    jsd_srs = pd.Series(index=align_df.columns)
    if keep_test_spec:
        jsd_df = align_df.copy()
    else:
        jsd_df = align_df.drop(test_idx,axis=0)
    jsd_nd = jsd_df.values
    weights = jsd_calc.calculate_sequence_weights(jsd_nd)
    for i,col in enumerate(jsd_nd.T):
        jsd = jsd_calc.JSD(col,blosum62_bg,weights,aas,use_gap_penalty=use_gap_penalty)
        jsd_srs[i+1] = jsd
    jsd_zscores = calc_z_scores(jsd_srs)
    return jsd_srs, jsd_zscores

def filter_score_srs(score_srs, uniques):
    """Filters score series down to positions in uniques.

    :param score_srs: pandas Series of score entries, indexed on alignment positions
    :param uniques: DataFrame corresponding to columns from alignment DataFrame with unique residue substitutions
    :return: unique_scores: pandas Series corresponding to entries in score_srs corresponding to alignment positions
    in uniques
    """
    unique_positions = uniques.columns
    unique_scores = score_srs.loc[unique_positions]
    return unique_scores

def variant_counts(col,test_spec_idx):
    """For col, return a dictionary containing information on gap_fraction,test_species/ outgroup variants, and the
    number of instances/ ratios of those variants

    :param (Series) col: Series corresponding to alignment dataframe column (including test species)
    :param test_spec_idx: index of col corresponding to test species record
    :return: metrics: dictionary containing test_variant, outgroup_variant, and gap_fraction
    """
    outgroup_variant = col.mode().iloc[0]
    gap_ratio = gap_fraction(col)
    test_variant = col[test_spec_idx][0]
    vc = col.value_counts()
    og_var_freq,test_var_freq = vc[outgroup_variant],vc[test_variant]
    og_var_ratio,test_var_ratio = og_var_freq/len(col), test_var_freq/len(col)
    vc_metrics = {'test_variant':test_variant,'outgroup_variant':outgroup_variant,'gap_fraction':gap_ratio,
               'outgroup_variant_count':og_var_freq,'test_variant_count':test_var_freq,
               'outgroup_variant_ratio':og_var_ratio,'test_variant_ratio':test_var_ratio}
    return vc_metrics

def test_outgroup_blosum(col,test_spec_idx,blos_df):
    """Calculates average BLOSUM62 score of the test species variant in col against all outgroup variants. NaN for gaps
    or other ambiguous IUPAC characters (ie B, U, X)

    :param col: pandas Series corresponding to column from align_df
    :param test_spec_idx: index corresponding to test species
    :param blos_df: DataFrame corresponding to blosum matrix; index and columns are amino acid characters
    :return: Average test variant vs outgroup variant blosum scores for col.
    """
    test_var = col[test_spec_idx][0]
    outgroup_col = col.drop(test_spec_idx)
    outgroup_col = outgroup_col[outgroup_col.isin(NON_GAP_AAS)]
    og_col_nd = outgroup_col.values
    if test_var in NON_GAP_AAS:
        blos_mean = np.mean(blos_df[test_var][og_col_nd])
    else:
        blos_mean = np.nan
    return blos_mean

def test_outgroup_blosum_series(align_df,test_spec_idx,blos_df):
    """Returns series of scores/z-scores for test vs outgroup blosum values for entire align_df.

    :param align_df: MSA DataFrame
    :param test_spec_idx: Index object corresponding to test_species. Remaining species in align_df considered outgroup
    :param blos_df: Blosum62 DataFrame
    :return: blos_srs: Series of test vs outgroup BLOSUM62 scores
    :return: blos_z: Z-scores calculated for above series
    """
    blos_srs = pd.Series(index=align_df.columns)
    for pos in align_df.columns:
        aln_col = align_df.loc[:,pos]
        blos_mean = test_outgroup_blosum(aln_col,test_spec_idx,blos_df)
        blos_srs[pos] = blos_mean
    blos_srs.name = "Test-Outgroup BLOSUM62"
    blos_z = calc_z_scores(blos_srs)
    return blos_srs, blos_z

def pairwise_outgroup_blosum(col,test_spec_idx,blos_df):
    """Returns average pairwise blosum score between all non-gap residues in outgroup of column

    :param col: pandas Series corresponding to column from alignment DataFrame
    :param test_spec_idx: index of test_species in col
    :param blos_df: BLOSUM DataFrame from gen_blos_df
    :return: Average pairwise BLOSUM score over all residues in outgroup against each other
    """
    outgroup_col = col.drop(test_spec_idx)
    og_col = outgroup_col.values
    pairwise_scores = []
    skip_chars = ['-', 'X']
    for i, first in enumerate(og_col):
        for j, second in enumerate(og_col):
            # if i < j and not (first in skip_chars) and not second in skip_chars:
            if i < j and (first in NON_GAP_AAS) and second in NON_GAP_AAS:
                pairwise_scores.append(blos_df[first][second])
            else:
                pass
    if pairwise_scores:
        blos_pw = np.mean(pairwise_scores)
    else:
        blos_pw = np.nan
    return blos_pw


#Global variables for module
aas, blosum62_bg, blos_df, sim_matrix = gen_blos_df()

NON_GAP_AAS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
UNAMBIG_CHARS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','-']