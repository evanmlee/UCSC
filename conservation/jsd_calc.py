#JSDcalc.py - Provides all functions for calculating JSD using Capra and Singh's 2007 method.
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
from IPython.display import display

def calculate_sequence_weights(msa):
    """Calculate the sequence weights using the Henikoff '94 method for the given multiple sequence alignment. Briefly,
    sequence weights are based on sequence diversity at each position in the MSA.
    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    https://compbio.cs.princeton.edu/conservation/
    :param msa: np.ndarray object containing characters from given multiple sequence alignment.
    :return: seq_weights: array of sequence weights for given msa.
    """
    from collections import Counter
    n = len(msa)
    alignlen = len(msa[0])
    seq_weights = np.zeros((n))
    for col in msa.T:
        freqs = Counter(col)
        del freqs['-']
        num_observed = len(freqs)
        d_ = [freqs[aa] * num_observed for aa in col]
        inverse = [1 / d if d > 0 else 0 for d in d_]
        seq_weights += inverse
    seq_weights /= alignlen
    return seq_weights

def weighted_fcpseudo(col, seq_weights, aas, pc_amount=0.0000001):
    """Generates Pc distribution for column col. Distribution is weighted using seq_weights
    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    https://compbio.cs.princeton.edu/conservation/
    :param col:
    :param seq_weights:
    :param pc_amount:
    :return:
    """
    if len(seq_weights) != len(col):
        return [1]*len(col)
    else:
        freqs = np.array([pc_amount]*len(aas))
        for i,aa in enumerate(aas):
            freqs[i] += sum([seq_weights[j] if aa == entry else 0 for j,entry in enumerate(col)])
        for i,entry in enumerate(freqs):
            freqs[i] /= (sum(seq_weights) + len(aas)*pc_amount)
        return freqs

def weighted_gap_penalty(col,weights):
    """Calculate the simple gap penalty multiplier for the column. If the sequences are weighted,
    the gaps, when penalized, are weighted accordingly.

    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    :param col:
    :param weights:
    :return:
    """
    gap_weights = np.array([w if c == "-" else 0 for c,w in zip(col,weights)])
    gap_sum = sum(gap_weights)
    return 1 - (gap_sum/sum(weights))
# print(weighted_gap_penalty(['A','A','-'],[0.3333,0.3333,0.3333]))
def relative_entropy(pC,bg,gap_penalty=1):
    """Calculate the relative entropy of the column distribution with a
    background distribution specified in bg_distr. This is similar to the
    approach proposed in Wang and Samudrala 06.
    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    https://compbio.cs.princeton.edu/conservation/
    :param pC: Column distribution of amino acids using fc_pseudocount
    :param bg: background distribution
    :param gap_penalty: If provided, will be used to adjust relative entropy. Calculated in JSD using weighted_gap_penalty
    """
    if len(bg) == 20 and len(pC) == 21:
        #Remove gap count
        pC = pC[:-1]
        pC = pC/(sum(pC))
    ratio = pC/bg
    log_ratio = np.log(ratio)/np.log(2)
    RE = sum([a*b for a,b in zip(pC,log_ratio)])
    RE*= gap_penalty
    return RE

def JSD(col,bg,weights,aas,jsd_lam=0.5,use_gap_penalty=True):
    """Jensen-Shannon Divergence, calculation based on below.
    John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics,
    Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270
    https://compbio.cs.princeton.edu/conservation/
    :param col: column distribution to be compared against bg
    :param bg: background distribution used in relative entropy calculation
    :param weights: sequence weights
    :param aas: amino acid characters, see SSanalysiscalc.py
    :param jsd_lam: 0.5 by default, lambda value to be used in jsd calculation (see paper)
    :return:
    """

    pC = weighted_fcpseudo(col,weights,aas)
    if len(bg) == 20:
        #Remove gap count
        pC = pC[:-1]
        pC = pC/(sum(pC))
    if use_gap_penalty:
        gap_penalty = weighted_gap_penalty(col,weights)
    else:
        gap_penalty = 1
    r_dist = np.array([jsd_lam*aa for aa in pC]) + np.array([(1-jsd_lam)*aa for aa in bg])

    RE_pCr = relative_entropy(pC,r_dist,gap_penalty)
    RE_qr = relative_entropy(bg,r_dist, gap_penalty)
    jsd = jsd_lam*RE_pCr + (1-jsd_lam)*RE_qr
    return jsd