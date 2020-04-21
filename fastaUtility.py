from Bio import Seq, SeqIO
import subprocess
import pandas as pd
import numpy as np
import os
import re
from collections import OrderedDict

from SSfasta import filter_fasta_infile

def profile_MSA(in1_fpath, in2_fpath, out_fpath):
    subprocess.run(args=['./muscle', '-profile', '-in1', in1_fpath,
                         '-in2', in2_fpath, '-out', out_fpath])

def fpath_list_record_generator(fpath_list):
    for fpath in fpath_list:
        fasta_seqs = SeqIO.parse(open(fpath), 'fasta')
        for fasta in fasta_seqs:
            yield fasta

def NCBI_ortholog_fpath_list(NCBI_parent_dir,NCBI_gid, tax_id_dict):
    """Given an NCBI data parent, checks for existing ortholog gene files, stored as [species_tax_id]/[NCBI_human_gid]

    :param NCBI_gid:
    :param tax_id_dict:
    :return: list of existing file paths matching that architecture
    """
    fpath_list = []
    for taxid in tax_id_dict:
        NCBI_fpath = "{0}/{1}/{2}.fasta".format(NCBI_parent_dir, taxid, NCBI_gid)
        if (os.path.exists(NCBI_fpath)):
            fpath_list.append(NCBI_fpath)
    return fpath_list

def MSA_fpath_list(fpath_list, combined_outpath, aligned_outpath):
    combined_generator = fpath_list_record_generator(fpath_list)
    SeqIO.write(combined_generator, combined_outpath, "fasta")
    subprocess.run(args=['./muscle', '-in', combined_outpath,
                         '-out', aligned_outpath])

def fasta_to_srs(fasta_path):
    fasta_seqs = SeqIO.parse(open(fasta_path), 'fasta')
    id_seq_map = OrderedDict()
    for fasta in fasta_seqs:
        record_id = fasta.id
        seq = str(fasta.seq)
        id_seq_map[record_id] = seq
    return pd.Series(name="seq", data=id_seq_map)

def construct_id_dm(seq_df, seq_fpath, align_outpath="tmp/iddm_align.fasta", ordered=False, aligned=False):
    """Constructs an np.ndarray corresponding to the identity distance matrix of records in seq_df. Different from
    SSfasta version, uses MUSCLE for alignment and allows for aligned parameter to prevent unnecessary MUSCLE calls.

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
    filtered_fpath = "tmp/iddm.fasta"
    filter_fasta_infile(seq_df.index, seq_fpath, outfile_path=filtered_fpath, ordered=ordered)
    if not aligned:
        # MUSCLE align sequences from filtered_fpath to align_outpath
        subprocess.run(args=['./muscle', '-in', filtered_fpath, \
                             '-out', align_outpath])
    else:
        #Use filtered file from filter_fasta_infile as input for identity distmat
        align_outpath = filtered_fpath
    align_srs = fasta_to_srs(align_outpath)
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

def min_dist_spec_record(distmat,dm_record_ids,spec_record_ids,accepted_record_ids, ref_df):
    """Calculates the average distance of every record containing spec_taxid against accepted records, then
    returns the row from ref_df corresponding to the record with lowest average distance.

    :param distmat: np.ndarray n x n distance matrix
    :param dm_record_ids: iterable of record_ids corresponding to rows of distmat
    :param spec_record_ids: list of records from specific species
    :param accepted_record_ids: record ids against which record distances will be calculated
    :param ref_df: Sequence dataframe
    :return: md_row: row from ref_df with minimum average distance to accepted_record_id records,
    :return min_dist: minimum average distance value
    """
    spec_records = [(i,id_) for i,id_ in enumerate(dm_record_ids) if id_ in spec_record_ids]
    spec_dm_idxs = [t[0] for t in spec_records]
    accepted_records = [(i,id_) for i,id_ in enumerate(dm_record_ids) if id_ in accepted_record_ids]
    accepted_dm_idxs = [t[0] for t in accepted_records]

    spec_dm = distmat[:,spec_dm_idxs]
    sub_dm = spec_dm[accepted_dm_idxs,:]
    if len(sub_dm) > 1:
        avg_dist = sub_dm.mean(axis=0)
    else:
        avg_dist = sub_dm[0]
    min_idx = np.argmin(avg_dist)
    min_dist_id = spec_records[min_idx][1]
    min_dist = avg_dist[min_idx]
    md_row = ref_df.loc[min_dist_id,:]
    return md_row, min_dist

###Fasta File into DataFrame objects; parses description fields into species_name, taxid; indexed by record_id###

def NCBI_fasta_df(fasta_inpath,tax_id_dict):
    #Dictionary of species name to tax_id
    spec_name_dict = {tax_id_dict[tid]: tid for tid in tax_id_dict}
    exact_spec_name_re = "({0})".format("|".join(spec_name_dict.keys()))
    fastas = SeqIO.parse(fasta_inpath,format='fasta')
    column_list = ["description","species_name","NCBI_taxid","length","sequence"]
    fasta_df = pd.DataFrame(columns=column_list)
    for fasta in fastas:
        desc_remaining = re.search("{0}(.*)".format(fasta.id),fasta.description).groups()[0]
        if desc_remaining:
            desc_remaining = desc_remaining.strip()
        else:
            desc_remaining = ''
        fasta_seq = str(fasta.seq)
        #fasta_length functional for both aligned and unaligned sequences
        fasta_length = len(fasta_seq) - fasta_seq.count("-")
        #matched_spec_name: Pattern match from NCBI sequence data - can contain third taxonomic subspecies word
        matched_spec_name = re.search("\[(\w+\s\w+(\s\w+)?)\]$",desc_remaining)
        if matched_spec_name:
            matched_spec_name = matched_spec_name.groups()[0]
        else:
            continue
        #spec_name: Exact string match to
        spec_name = re.search(exact_spec_name_re,matched_spec_name).groups()[0]
        taxid = spec_name_dict[spec_name]

        row_vals = [desc_remaining,spec_name,taxid,fasta_length,fasta_seq]
        row = dict(zip(column_list,row_vals))
        fasta_df.loc[fasta.id,:] = row
    fasta_df.index.name = "record_id"
    return fasta_df

def UCSC_fasta_df(fasta_inpath):
    fastas = SeqIO.parse(fasta_inpath, format='fasta')
    column_list = ["description", "species_name", "length","sequence"]
    fasta_df = pd.DataFrame(columns=column_list)
    for i,fasta in enumerate(fastas):
        desc_remaining = re.search("{0}(.*)".format(fasta.id), fasta.description).groups()[0].strip()
        #Remove UCSC specific stop character (Z at last position) 
        fasta_seq = str(fasta.seq)[:-1]
        #Remove gap characeters and UCSC stop character (Z) from align_length to get sequence length
        fasta_length = len(fasta_seq) - fasta_seq.count("Z") - fasta_seq.count("-")
        spec_name = re.search("_(\w+)$",fasta.id).groups()[0]
        row_vals = [desc_remaining, spec_name, fasta_length, fasta_seq]
        row = dict(zip(column_list, row_vals))
        fasta_df.loc[fasta.id, :] = row
    fasta_df.index.name = "record_id"
    return fasta_df

def select_NCBI_records(dir_vars,tax_id_dict,UCSC_tid,NCBI_gid,selection="identity"):
    from IPython.display import display
    UCSC_parent,allNCBI_parent,bestNCBI_parent = [dir_vars[k] for k in ["UCSC_raw_parent","allNCBI_parent",
                                                                        "bestNCBI_parent"]]

    raw_UCSC_fpath = "{0}/{1}.fasta".format(UCSC_parent,UCSC_tid)
    NCBI_all_aln_fpath = "{0}/NCBI_alignments/{1}.fasta".format(allNCBI_parent,NCBI_gid)
    combined_all_aln_fpath = "{0}/aligned/{1}.fasta".format(allNCBI_parent,UCSC_tid)

    UCSC_df = UCSC_fasta_df(raw_UCSC_fpath)
    allNCBI_df = NCBI_fasta_df(NCBI_all_aln_fpath,tax_id_dict)

    CLOSEST_EVO = {'9999':['mm10','rn6','speTri2'],'29073':['ailMel1','canFam3'],
                   '10181':['hetGla2','mm10','rn6'], '9994':['mm10','rn6','speTri2']}

    allseq_df = UCSC_df.append(allNCBI_df,sort=False)
    allseq_df.loc[allNCBI_df.index,"NCBI_taxid"] = allNCBI_df["NCBI_taxid"]
    # display(allseq_df["sequence"])
    id_dm, align_srs = construct_id_dm(allseq_df,combined_all_aln_fpath,aligned=True)
    bestNCBI_ids = []
    bestNCBI_dict = {}
    for taxid in tax_id_dict:
        NCBI_tid_df = allNCBI_df.loc[allNCBI_df["NCBI_taxid"]==taxid,:]
        #Skip if no taxid records in allNCBI_df (ie no records for that species)
        if not NCBI_tid_df.empty:
            spec_record_ids = NCBI_tid_df.index
            closest_evos = CLOSEST_EVO[taxid]
            ce_pat = '|'.join(closest_evos)
            id_calc_records = UCSC_df.index[UCSC_df.index.str.contains(ce_pat)]
            #Check if records in id_calc_records are all empty strings
            if sum([len(seq) for seq in UCSC_df.loc[id_calc_records,"sequence"]]) == 0:
                #TODO: Figure out filtering if no UCSC records for comparison (likely only for polar bear NCBI)
                continue
            md_row, min_dist = min_dist_spec_record(id_dm,align_srs.index,spec_record_ids,id_calc_records,allseq_df)
            bestNCBI_ids.append(md_row.name)
            bestNCBI_dict[taxid] = md_row

    bestNCBI_df = allNCBI_df.loc[bestNCBI_ids]
    # display(bestNCBI_df)
    bestNCBI_aln_fpath = "{0}/NCBI_alignments/{1}.fasta".format(bestNCBI_parent,NCBI_gid)
    bestNCBI_all_fpath = "{0}/combined/{1}.fasta".format(bestNCBI_parent,UCSC_tid)
    if True:#not os.path.exists(bestNCBI_aln_fpath): (os.path check now in outer function in UCSC_filter)
        filter_fasta_infile(bestNCBI_ids,NCBI_all_aln_fpath,bestNCBI_aln_fpath)
        final_align_ids = list(UCSC_df.index)
        final_align_ids.extend(bestNCBI_ids)
        # print(final_align_ids)
        filter_fasta_infile(final_align_ids,combined_all_aln_fpath,bestNCBI_all_fpath)
    #UCSC_euth: Records which taxonomically belong to euarchontoglires or laurasiatheria
    UCSC_euth_df = UCSC_df.iloc[:51, :]
    UCSC_euth_mean_len, UCSC_euth_median_len = UCSC_euth_df["length"].mean(), UCSC_euth_df["length"].median()
    allNCBI_mean_len, allNCBI_median_len = allNCBI_df["length"].mean(),allNCBI_df["length"].median()
    bestNCBI_mean_len, bestNCBI_median_len = bestNCBI_df["length"].mean(), bestNCBI_df["length"].median()
    length_labels = ["euth_mean","euth_median","allncbi_mean","allncbi_median","bestncbi_mean","bestncbi_median"]
    length_vals = [UCSC_euth_mean_len,UCSC_euth_median_len,
                   allNCBI_mean_len,allNCBI_median_len,
                   bestNCBI_mean_len,bestNCBI_median_len]
    length_metrics = dict(zip(length_labels,length_vals))
    return length_metrics
