from Bio import Seq, SeqIO, SeqRecord
from Bio.Alphabet import IUPAC
import subprocess
import pandas as pd
import numpy as np
import os
import re
from collections import OrderedDict
import warnings

###Record Filtering###

def ordered_record_generator(fpath,ordered_ids):
    """Generates Seq records from fpath, limited to ordered_ids if provided. Generator order is order in ordered_ids.
    Internally tracks yielded record ids to eliminate duplicates (keeps first duplicate only).
    :param fpath: Fasta file path
    :param ordered_ids: array-like containing ordered record ids
    :return:
    """
    yielded = set()
    for id in ordered_ids:
        with open(fpath) as fasta_f:
            fastas = SeqIO.parse(fasta_f,'fasta')
            for fasta in fastas:
                if fasta.id == id and fasta.id not in yielded:
                    yielded.add(fasta.id)
                    yield fasta
                    break

def record_generator(fpath,ids=[]):
    """Generates Seq records from fpath, limited to ids if provided. Ordered as in fpath.
    Internally tracks yielded record ids to eliminate duplicates (keeps first duplicate only).
    :param fpath: Fasta file path
    :param ids: If provided, only yield records corresponding to record ids provided. Default value causes all records
    to be present in generator
    :return: generator object
    """
    yielded = set()
    with open(fpath) as fasta_f:
        fastas = SeqIO.parse(fasta_f, 'fasta')
        for fasta in fastas:
            if ((len(ids)>0 and fasta.id in ids) or \
                    (len(ids)==0)) and fasta.id not in yielded:
                yielded.add(fasta.id)
                yield fasta

def UCSC_subtax_generator(UCSC_fpath,clade='all'):
    low,high = UCSC_clade_positions(clade)
    idx = 0
    for fasta in SeqIO.parse(UCSC_fpath,'fasta'):
        if idx >= low and idx < high:
            yield fasta
        idx += 1

def UCSC_NCBI_generator(UCSC_fpath,NCBI_fpath,ucsc_subset=[],ncbi_subset=[],ordered=False):
    """Generator object that yields all Bio.Seq objects in ODB_fpath and then NCBI_fpath, filtered down with optional
    parameters odb_subset and ncbi_subset.

    :param ODB_fpath: file path to OrthoDB fasta
    :param NCBI_fpath: file path to NCBI fasta
    :param odb_subset: If provided, only ODB records present in odb_subset will be provided by generator. Else all
    records will be provided
    :param ncbi_subset: If provided, only NCBI records present in ncbi_subset will be provided by generator. Else all
    records will be provided
    :param ordered: If true, odb_subset and ncbi_subset must be provided. Causes records to be returned in order they
    are present in odb_subset and ncbi_subset.
    :return: Generator object which provides Bio.Seq objects filtered and ordered as described above.
    """

    if ordered:
        if len(ucsc_subset) == 0 or len(ncbi_subset)==0:
            raise ValueError("Both odb_subset and ncbi_subset must be provided if ordered is True.")
        else:
            odb_generator = ordered_record_generator(UCSC_fpath,ucsc_subset)
            ncbi_generator = ordered_record_generator(NCBI_fpath,ncbi_subset)
    else:
        odb_generator = record_generator(UCSC_fpath, ucsc_subset)
        ncbi_generator = record_generator(NCBI_fpath, ncbi_subset)

    for fasta in odb_generator:
        yield fasta
    for fasta in ncbi_generator:
        yield fasta

def filtered_generator_wrapper(filtered_ids, infile_path, ordered, missing_warning=False):
    """Wrapper function which generates appropriate generatoe object based on ordered. Issues warnings for missing
    values from filtered_ids if not present in infile_path if missing_warning is True.

    :param filtered_ids: See filter_fasta_infile
    :param infile_path: See filter_fasta_infile
    :param ordered: if True, generator returns records in order of appearance in filtered_ids; False -> order from
    record order in infile.
    :param missing_warning: boolean on whether to issue warnings for missing filtered_id values from infile
    :return:
    """
    if missing_warning:
        with open(infile_path) as infile_f:
            infile_ids = [fasta.id for fasta in SeqIO.parse(infile_f, 'fasta')]
            for record_id in filtered_ids:
                if record_id not in infile_ids:
                    msg = "Infile {0} is missing record id {1} from filtered_ids".format(infile_path, record_id)
                    warnings.warn(msg)
    if ordered:
        # filtered = _ordered_filtered_generator(filtered_ids, infile_path)
        filtered = ordered_record_generator(infile_path,filtered_ids)
    else:
        filtered = record_generator(infile_path, filtered_ids)
        # filtered = _filtered_generator(filtered_ids, infile_path)
    return filtered

def filter_fasta_infile(filtered_ids, infile_path, outfile_path=None, ordered=False):
    """Filters fasta_infile records from infile_path and returns a series of fasta records if id in filtered_ids

    :param filtered_ids: Iterable containing all id values for which records will be saved
    :param infile_path: Fasta file containing unfiltered set of records. If missing records from filtered_ids,
    will issue warning
    :param outfile_path: Optional filepath parameter. If provided, filtered results will be written to file path
    :param ordered: boolean, if true, sequences will be returned/ written in order of filtered_ids
             if false, uses sequence order of sequences in infile_path
    :return: Series of fasta sequences (index is fasta.id) for which id is present in filtered_ids
    """
    if outfile_path:
        filtered = filtered_generator_wrapper(filtered_ids,infile_path,ordered)
        SeqIO.write(filtered, outfile_path, "fasta")
    with warnings.catch_warnings(record=True) as w:
        filtered = filtered_generator_wrapper(filtered_ids, infile_path, ordered, True)
    filtered_srs = pd.Series(index=filtered_ids)
    for fasta in filtered:
        filtered_srs[fasta.id] = str(fasta.seq)
    return filtered_srs

###Profile alignment and NCBI record compilation functions###


def profile_MSA(in1_fpath, in2_fpath, out_fpath):
    """Generates profile-profile alignment of records from in1_fpath and in2_fpath and writes to out_fpath
    :return: N/A
    """
    subprocess.run(args=['./muscle', '-profile', '-quiet','-in1', in1_fpath,
                         '-in2', in2_fpath, '-out', out_fpath])

def MSA(in_fpath,out_fpath,aln_program='muscle'):
    """Generates a multiple sequence alignment at out_fpath using the specified program."""
    if aln_program == 'muscle':
        subprocess.run(args=['./muscle','-quiet','-in', in_fpath, '-out', out_fpath])
    elif aln_program == 'kalign':
        args = ['kalign']
        subprocess.run(args=args, stdin=in_fpath, stdout=out_fpath, stderr=subprocess.PIPE, text=True)
    else:
        raise ValueError("Parameter 'aln_program' must be 'muscle' or 'kalign'.")

def fpath_list_record_generator(fpath_list):
    """Creates a generator object which yiellds all the sequences contained in the files in fpath_list.
    """
    for fpath in fpath_list:
        with open(fpath) as seq_f:
            fasta_seqs = SeqIO.parse(seq_f, 'fasta')
            for fasta in fasta_seqs:
                yield fasta

def NCBI_ortholog_fpath_list(NCBI_parent_dir,NCBI_gid,taxid_dict=None):
    """Given an NCBI data parent, checks for existing ortholog gene files, stored as [species_tax_id]/[NCBI_human_gid]

    :param NCBI_parent_dir: Parent orthologs sequence data directory path
    :param NCBI_gid: Human NCBI Gene ID to check for which to check corresponding ortholog files
    :param taxid_dict: If provided, will only use NCBI_orthologs represented in provided dict. Otherwise, will use
    taxids represented in dict read from ncbi_analysis_tax.txt
    :return: list of existing file paths (ie list of species-specific sequence data for which ortholog data
    was available and successfully downloaded from NCBI.
    """
    fpath_list = []
    if not taxid_dict:
        taxid_dict = table_taxid_dict
    for taxid in taxid_dict:
        NCBI_fpath = "{0}/{1}/{2}.fasta".format(NCBI_parent_dir, taxid, NCBI_gid)
        if (os.path.exists(NCBI_fpath)):
            fpath_list.append(NCBI_fpath)
    return fpath_list

def MSA_fpath_list(fpath_list, combined_outpath, aligned_outpath,aln_program='muscle'):
    """Generates both combined unaligned and aligned files for NCBI records in fpath_list

    :param fpath_list: List of file paths to Ortholog gene files
    :param combined_outpath: Outpath to write combined (unaligned) NCBI records
    :param aligned_outpath: Outpath to write combined aligned NCBI records
    :param aln_program:
    :return: N/A. Generates appropriate files at combined_outpath and aligned_outpath
    """
    combined_generator = fpath_list_record_generator(fpath_list)
    SeqIO.write(combined_generator, combined_outpath, "fasta")
    MSA(combined_outpath,aligned_outpath,aln_program)

def construct_id_dm(record_ids, seq_fpath, align_outpath="tmp/iddm_align.fasta", ordered=False,filtered=False,aligned=False):
    """Constructs an np.ndarray corresponding to the identity distance matrix of records in seq_df. Different from
    SSfasta version, uses MUSCLE for alignment and allows for aligned parameter to prevent unnecessary MUSCLE calls.

    :param record_ids: Record IDs for which identity distance matrix will be computed
    :param seq_fpath:  Path of fasta file containing at least all of the records in seq_df. Can contain more records
    than are in seq_df - a temporary file containing only the records in seq_df.index will be generated (filtered_fpath)
    :param align_outpath: Optional filepath. If provided, the resulting alignment will be stored there. Otherwise,
    written to a temporary file (tmp/iddm_align.fasta)
    :param ordered: boolean. True: distance matrix rows will be ordered by the order of records in seq_df.index;
    False: distance matrix rows will be ordered by the order of records in seq_fpath
    :param: (boolean) filtered: If True, unfiltered record set from seq_fpath will be used for alignment and dm calc.
    :param: (boolean) aligned: If True, filtered records from seq_fpath must be aligned and ready for use in dm calculation.
    :return: id_dm: np.ndarray of identity distance matrix calculated by AlignIO
    :return: align_srs: pandas Series object containing aligned sequences
    """
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio import AlignIO

    # Filter records in seq_fpath to new fasta only containing records in seq_df.index
    if not filtered:
        filtered_fpath = "tmp/iddm.fasta"
        filter_fasta_infile(record_ids, seq_fpath, outfile_path=filtered_fpath, ordered=ordered)
    else:
        filtered_fpath = seq_fpath
    if not aligned:
        # MUSCLE align sequences from filtered_fpath to align_outpath
        MSA(filtered_fpath,align_outpath,aln_program='muscle')
    else:
        #Use filtered file from filter_fasta_infile as input for identity distmat
        align_outpath = filtered_fpath
    align_srs = fasta_to_srs(align_outpath)
    with open(align_outpath) as align_f:
        aln = AlignIO.read(align_outpath, 'fasta')
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

###Sequence File/Series/DataFrame Manipulations###

def srs_to_fasta(seq_srs, outfile_path):
    # Write records in seq_srs to outfile_path in fasta format
    #Deprecated; use record generator functions to write sequence data (maintains description/ other non-sequence info)
    def record_generator(seq_srs):
        for idx, seq in seq_srs.iteritems():
            record = SeqRecord(Seq(seq, IUPAC.protein), id=idx)
            yield record

    records = record_generator(seq_srs)
    SeqIO.write(records, outfile_path, "fasta")

def fasta_to_srs(fasta_path,ucsc_flag=False):
    #Creates series mapping record id to sequence from fasta_path
    with open(fasta_path) as fasta_f:
        fasta_seqs = SeqIO.parse(fasta_f, 'fasta')
        id_seq_map = OrderedDict()
        for fasta in fasta_seqs:
            record_id = fasta.id
            seq = str(fasta.seq)
            if ucsc_flag:
                seq = seq.replace('Z','-')
            id_seq_map[record_id] = seq
        return pd.Series(name="seq", data=id_seq_map)

def align_srs_to_df(align_srs):
    """Returns DataFrame object from series of aligned sequences; columns are 1-indexed positions
    Values are characters in alignment, indexed on record_ids"""
    seq_len = len(align_srs.iloc[0])
    align_df = pd.DataFrame(index=align_srs.index, columns=range(seq_len))
    for idx, seq in align_srs.iteritems():
        align_df.loc[idx, :] = list(seq)
    align_df.columns += 1
    return align_df

def align_fasta_to_df(fasta_path,ucsc_flag=False):
    align_srs = fasta_to_srs(fasta_path,ucsc_flag=ucsc_flag)
    align_df = align_srs_to_df(align_srs)
    return align_df

def seq_srs_to_align_df(seq_srs, align_in_fpath, align_out_fpath):
    """Transform seq_srs (pandas Series containing sequence texts) to a DataFrame for which each column
    is an alignment position and column. Writes input fasta and output fastas for alignment to align_in_fpath
    and align_out_fpath respectively. Also returns average (non-diagonal) identity distances"""
    srs_to_fasta(seq_srs, align_in_fpath)
    id_dm, align_srs = construct_id_dm(seq_srs.index, align_in_fpath, align_out_fpath)
    align_df = align_srs_to_df(align_srs)
    return align_df

def align_srs_to_seq_srs(align_srs, outfile_path=None):
    """Return new Series (same index) of sequences with gap characters dropped
    If outfile_path is provided, write un-aligned record seqs to new fasta file"""
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


###Fasta File into record DataFrame objects; parses description fields into species_name, taxid; indexed by record_id###

def NCBI_fasta_df(fasta_inpath,taxid_dict=None):
    #Dictionary of species name to tax_id
    if not taxid_dict:
        taxid_dict = table_taxid_dict
    tids, spec_names = list(taxid_dict.keys()),[taxid_dict[k] for k in taxid_dict]
    exact_spec_name_re = "({0})".format("|".join(spec_names))
    spec_name_dict = dict(zip(spec_names,tids))
    fastas = SeqIO.parse(fasta_inpath,format='fasta')
    column_list = ["description","species_name","NCBI_taxid","length","sequence"]
    fasta_df = pd.DataFrame(columns=column_list)
    for fasta in fastas:
        ncbi_row = parse_NCBI_row(fasta,exact_spec_name_re,spec_name_dict)
        fasta_df.loc[fasta.id,:] = ncbi_row
    fasta_df.index.name = "record_id"
    return fasta_df

def parse_NCBI_row(fasta,exact_spec_name_re,spec_name_dict):
    """Converts Bio.SeqRecord object into DataFrame row populated with NCBI record data.
    :param fasta: Bio.SeqRecord object corresponding to fasta record
    :param exact_spec_name_re: regular expression pattern matching all species names in spec_name_dict
    :return: ncbi_row: pandas Series object with populated record information fields
    """
    column_list = ["description", "species_name", "NCBI_taxid", "length", "sequence"]
    desc_remaining = re.search("{0}(.*)".format(fasta.id), fasta.description).groups()[0]
    if desc_remaining:
        desc_remaining = desc_remaining.strip()
    else:
        desc_remaining = ""
    fasta_seq = str(fasta.seq)
    # fasta_length functional for both aligned and unaligned sequences
    fasta_length = len(fasta_seq) - fasta_seq.count("-")
    # matched_spec_name: Pattern match from NCBI sequence data - can contain third taxonomic subspecies word
    spec_name_match = re.search(r"\[(\w+\s\w+(\s\w+)?)\]$", desc_remaining)
    if spec_name_match:
        matched_name = spec_name_match.groups()[0]
        # spec_name: Exact string match to
        exact_name_match = re.search(exact_spec_name_re, matched_name)
        if not exact_name_match:
            spec_name, taxid = "",""
        else:
            spec_name = exact_name_match.groups()[0]
            taxid = spec_name_dict[spec_name]
    else:
        spec_name,taxid = "",""
    row_vals = [desc_remaining, spec_name, taxid, fasta_length, fasta_seq]
    row = dict(zip(column_list, row_vals))
    return pd.Series(data=row,name=fasta.id)

def UCSC_fasta_df(fasta_inpath):
    fastas = SeqIO.parse(fasta_inpath, format='fasta')
    column_list = ["description", "species_name", "NCBI_taxid", "length","sequence"]
    fasta_df = pd.DataFrame(columns=column_list)
    for i,fasta in enumerate(fastas):
        ucsc_row = parse_UCSC_row(fasta)
        fasta_df.loc[fasta.id, :] = ucsc_row
    fasta_df.index.name = "record_id"
    return fasta_df

def UCSC_clade_positions(clade):
    """Maps clade string identifier to low,high 0-indexed positions in UCSC tax_table or alignments for use in slices"""
    clade = clade.lower()
    clades = ['all','mammalia', 'boreoeutheria', 'euarchontoglires', 'primates']
    pos_tups = [(0,100),(0, 62), (0, 51), (0, 26), (0, 12)]
    idx_dict = dict(zip(clades, pos_tups))
    if clade not in clades:
        raise ValueError("Please specify accepted clade value (mammalia, boreoeutheria, euarchontoglires, primates.")
    else:
        low, high = idx_dict[clade]
    return low,high

def partition_UCSC_by_clade(fasta_df,clade):
    """
    :param fasta_df: Records DataFrame to filter data from
    :param (str) clade:
    :return:
    """
    low,high = UCSC_clade_positions(clade)
    incl_df,rest_df = partition_UCSC_by_position(fasta_df,low,high)
    return incl_df, rest_df

def partition_UCSC_by_position(fasta_df,low=0,high=100):
    """Partition fasta_df based on int_idx positions in the UCSC_tax_table. Returns incl and rest record DataFrames."""
    #Filter DataFrame to only UCSC formatted records
    fasta_df = fasta_df.loc[fasta_df.index.str.contains("ENST"),:]
    ucsc_tax_table_incl = ucsc_tax_table.iloc[low:high, :]
    ucsc_tax_table_rest = ucsc_tax_table.iloc[high:, :]
    incl_tids, rest_tids = ucsc_tax_table_incl['NCBI_taxid'], ucsc_tax_table_rest['NCBI_taxid']
    incl_df = fasta_df.loc[fasta_df['NCBI_taxid'].isin(incl_tids), :]
    rest_df = fasta_df.loc[fasta_df['NCBI_taxid'].isin(rest_tids), :]
    return incl_df, rest_df


def parse_UCSC_row(fasta):
    column_list = ["description", "species_name", "NCBI_taxid","length", "sequence"]
    desc_remaining = re.search("{0}(.*)".format(fasta.id), fasta.description).groups()[0].strip()
    # Remove UCSC specific stop character (Z at last position)
    fasta_seq = str(fasta.seq)[:-1]
    # Remove gap characeters and UCSC stop character (Z) from align_length to get sequence length
    fasta_length = len(fasta_seq) - fasta_seq.count("Z") - fasta_seq.count("-")
    fasta_handle_match = re.search("_(\w+)$", fasta.id)
    if fasta_handle_match:
        fasta_handle = fasta_handle_match.groups()[0]
    else:
        raise ValueError("Couldn't extract genome identifier from fasta.id: {0}".format(fasta.id))
    tax_table_row = ucsc_tax_table.loc[ucsc_tax_table['fasta_handle']==fasta_handle,:]

    spec_name, ncbi_taxid = tax_table_row['tax_name'].iloc[0], tax_table_row['NCBI_taxid'].iloc[0]
    row_vals = [desc_remaining, spec_name, ncbi_taxid,fasta_length, fasta_seq]
    row = dict(zip(column_list, row_vals))
    return pd.Series(data=row,name=fasta.id)

def load_UCSC_NCBI_df(combined_fpath,ncbi_taxid_dict={},UCSC_subset=[],NCBI_subset=[]):
    if not ncbi_taxid_dict:
        ncbi_taxid_dict = table_taxid_dict
    tids, spec_names = list(ncbi_taxid_dict.keys()),[ncbi_taxid_dict[k] for k in ncbi_taxid_dict]
    exact_spec_name_re = "({0})".format("|".join(spec_names))
    spec_name_dict = dict(zip(spec_names, tids))
    column_list = ["description", "species_name", "NCBI_taxid", "length", "sequence"]
    fasta_df = pd.DataFrame(columns=column_list)

    fastas = SeqIO.parse(combined_fpath, format='fasta')
    fasta_ids = [fasta.id for fasta in fastas]

    #If subset variables are not provided, default UCSC/ NCBI record partitioning
    if len(UCSC_subset) == 0:
        UCSC_subset = [fid for fid in fasta_ids if fid.__contains__('ENST')]
    if len(NCBI_subset) == 0:
        NCBI_subset =  [fid for fid in fasta_ids if not fid.__contains__('ENST')]
    fastas = SeqIO.parse(combined_fpath, format='fasta')
    for fasta in fastas:
        if fasta.id in UCSC_subset:
            ucsc_row = parse_UCSC_row(fasta)
            fasta_df.loc[fasta.id, :] = ucsc_row
        elif fasta.id in NCBI_subset:
            ncbi_row = parse_NCBI_row(fasta,exact_spec_name_re,spec_name_dict)
            fasta_df.loc[fasta.id, :] = ncbi_row
    fasta_df.index.name = "record_id"
    return fasta_df

def load_UCSC_tax_table(tsv_fpath="config/UCSC_tax.tsv"):
    tax_table = pd.read_csv(tsv_fpath,sep='\t',index_col='UCSC_spec_idx')
    #Extract genome identifiers used in fasta files from assembly version info
    fasta_handles = tax_table['assembly_version'].str.extract(r'/(\w+)$')
    tax_table['fasta_handle'] = fasta_handles
    return tax_table

def load_NCBI_tax_table(tsv_fpath="config/NCBI_analysis_tax.txt"):
    tax_table = pd.read_csv(tsv_fpath,sep='\t',index_col='spec_id')
    table_taxid_dict = {row['NCBI_taxid']: row['spec_name'] for i, row in tax_table.iterrows()}
    return tax_table, table_taxid_dict

ucsc_tax_table = load_UCSC_tax_table()
ncbi_tax_table, table_taxid_dict = load_NCBI_tax_table()
