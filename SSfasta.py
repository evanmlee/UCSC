from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from collections import OrderedDict
import pandas as pd
import numpy as np
import subprocess
import SSfasta
# Fasta file reading functions:
# filter_fasta_infile reads input files and outputs all records corresponding to filtered_ids to a new file
# Remaining functions provide conversions between fasta files, pandas Series, and pandas dataframes
# having alignment positions as columns

def filter_fasta_infile(filtered_ids, infile_path, outfile_path=None, ordered=False):
    # If outfile_path is provided, write filtered fasta to outfile_path
    """Generates new fasta file to outfile_path using the subset of sequences in infile_path
    which have ids in filtered_ids
    ordered: if true, sequences will be returned/ written in order of filtered_ids
             if false, uses sequence order of sequences in infile_path
    """

    def filtered_generator(filtered_ids, infile_path):
        fasta_seqs = SeqIO.parse(open(infile_path), 'fasta')
        for fasta in fasta_seqs:
            if fasta.id in filtered_ids:
                yield fasta

    def ordered_filtered_generator(filtered_ids, infile_path):
        for id_ in filtered_ids:
            fasta_seqs = SeqIO.parse(open(infile_path), 'fasta')
            for fasta in fasta_seqs:
                if fasta.id == id_:
                    yield fasta
                    break

    if outfile_path:
        if ordered:
            filtered = ordered_filtered_generator(filtered_ids, infile_path)
        else:
            filtered = filtered_generator(filtered_ids, infile_path)
        SeqIO.write(filtered, outfile_path, "fasta")
    if ordered:
        filtered = ordered_filtered_generator(filtered_ids, infile_path)
    else:
        filtered = filtered_generator(filtered_ids, infile_path)
    filtered_srs = pd.Series(index=filtered_ids)
    for fasta in filtered:
        filtered_srs[fasta.id] = str(fasta.seq)
    return filtered_srs


def srs_to_fasta(seq_srs, outfile_path):
    # Write records in seq_srs to outfile_path in fasta format
    def record_generator(seq_srs):
        for idx, seq in seq_srs.iteritems():
            record = SeqRecord(Seq(seq, IUPAC.protein), id=idx)
            yield record

    records = record_generator(seq_srs)
    SeqIO.write(records, outfile_path, "fasta")


def fasta_to_srs(fasta_path):
    fasta_seqs = SeqIO.parse(open(fasta_path), 'fasta')
    id_seq_map = OrderedDict()
    for fasta in fasta_seqs:
        record_id = fasta.id
        seq = str(fasta.seq)
        id_seq_map[record_id] = seq
    return pd.Series(name="seq", data=id_seq_map)


def align_srs_to_df(align_srs):
    # Returns DataFrame object from series of aligned sequences; columns are 1-indexed positions
    # Values are characters in alignment, index is ODB sequence IDs
    n_seq = len(align_srs)
    #     display(align_srs)
    #     display(align_srs.iloc[0])
    seq_len = len(align_srs.iloc[0])
    align_df = pd.DataFrame(index=align_srs.index, columns=range(seq_len))
    for idx, seq in align_srs.iteritems():
        align_df.loc[idx, :] = list(seq)
    align_df.columns += 1
    return align_df


def seq_srs_to_align_df(seq_srs, align_in_fpath, align_out_fpath):
    """Transform seq_srs (pandas Series containing sequence texts) to a DataFrame for which each column
    is an alignment position and column. Writes input fasta and output fastas for alignment to align_in_fpath
    and align_out_fpath respectively. Also returns average (non-diagonal) identity distances"""
    srs_to_fasta(seq_srs, align_in_fpath)
    n, ordered_ids, id_dm, align_srs = construct_id_dm(seq_srs, align_in_fpath, align_out_fpath)
    align_df = align_srs_to_df(align_srs)
    # dist_srs = avg_dist_srs(align_srs.index, id_dm)
    return align_df#, dist_srs


def align_srs_to_seq_srs(align_srs, outfile_path=None):
    # Return new Series (same index) of sequences with gap characters dropped
    # If outfile_path is provided, write un-aligned record seqs to new fasta file
    seq_srs = pd.Series(index=align_srs.index)
    for idx, align_seq in align_srs.iteritems():
        seq = align_seq.replace("-", "")
        seq_srs[idx] = seq
    if outfile_path:
        srs_to_fasta(seq_srs, outfile_path)
    return seq


def align_df_to_srs(align_df):
    # Returns series of aligned sequences from array of aligned positions
    align_srs = pd.Series(index=align_df.index)
    for idx, record in align_df.iterrows():
        #       #seq is a string joining all characters with no delimiter (i.e. the original aligned sequence with gaps)
        seq = ''.join(record.values)
        align_srs[idx] = seq
    return align_srs

### Distance Matrix Functions
#construct_id_dm makes an np.ndarray for the identity distance matrix of sequences for which OrthoDB id is
# in the index of seq_df; distance matrix rows will be ordered to the order in seq_fpath if ordered isFalse
def construct_id_dm(seq_df, seq_fpath, align_outpath="tmp/iddm_align.fasta", ordered=False):
    """Constructs an np.ndarray corresponding to the identity distance matrix of records in seq_df

    :param seq_df: DataFrame of OrthoDB/ NCBI sequence records; should only contain records for which identity
    distance matrix will be computed
    :param seq_fpath:  Path of fasta file containing at least all of the records in seq_df. Can contain more records
    than are in seq_df - a temporary file containing only the records in seq_df.index will be generated (filtered_fpath)
    :param align_outpath: Optional filepath. If provided, the resulting alignment will be stored there. Otherwise,
    written to a temporary file (tmp/iddm_align.fasta)
    :param ordered: boolean. True: distance matrix rows will be ordered by the order of records in seq_df.index;
    False: distance matrix rows will be ordered by the order of records in seq_fpath
    :return: id_dm: np.ndarray of identity distance matrix calculated by AlignIO
    :return: align_srs: pandas Series object containing aligned sequences
    """
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio import AlignIO
    # Filter records in seq_fpath to new fasta only containing records in seq_df.index
    # filtered_outpath = "tmp/iddm.fasta"
    filtered_fpath = "tmp/alias_matches.fasta"
    SSfasta.filter_fasta_infile(seq_df.index, seq_fpath, outfile_path=filtered_fpath, ordered=ordered)
    # KAlign sequences in filtered_outpath, write to align_outpath
    # n, ordered_ids, ka_dm, align_outfile = load_ka_distmat(filtered_outpath, align_outfile=align_outpath)
    subprocess.run(args=["kalign", '-i', filtered_fpath, "-o", align_outpath, "-f", "fasta"])
    align_srs = SSfasta.fasta_to_srs(align_outpath)
    aln = AlignIO.read(open(align_outpath), 'fasta')
    calculator = DistanceCalculator('identity')
    id_dm_obj = calculator.get_distance(aln)
    # Convert AlignIO object to np.ndarray
    for i, r in enumerate(id_dm_obj):
        if i == 0:
            id_dm = np.array(r)
        else:
            id_dm = np.vstack((id_dm, r))
    return id_dm, align_srs

def avg_dist_srs(index,distmat):
    #index is a pandas Index object with entries corresponding to the distmat (i.e. lengths should be equal)
    #Calculate mean of non-self record distances (diagonal distances generally force-set to 0, so
    #sum functions as intended)
    n = len(distmat)
    avg_dists = np.sum(distmat, axis=1)/(n-1)
    dist_srs = pd.Series(data=avg_dists,index=index,name="dist")
    return dist_srs