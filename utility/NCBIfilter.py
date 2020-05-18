import pandas as pd
import numpy as np
import os
from Bio import Seq, SeqIO
from collections import OrderedDict

from utility import directoryUtility, UCSCerrors
from IPython.display import display
from utility.fastaUtility import UCSC_fasta_df, NCBI_fasta_df, filter_fasta_infile, min_dist_spec_record, \
    construct_id_dm,ucsc_tax_table
from utility import fastaUtility as fautil
from query import orthologUtility as orutil

def load_transcript_table(transcripts_fpath):
    """Loads generic transcript ID indexed table into DataFrame."""
    transcript_table = pd.read_csv(transcripts_fpath,sep='\t',index_col="UCSC_transcript_id")
    return transcript_table

def read_evo_relative_data():
    """Reads closest and secondary evolutionary relative data specified in config/config.txt, returns dictionaries
    mapping NCBI anlaysis species Taxonomy IDs to either NCBI Tax IDs or UCSC record handles of corresponding
    evolutionary relative species in UCSC record data.
    :return: closest_evo_taxids,closest_evo_handles,second_evo_taxids,second_evo_handles: Map NCBI species taxonomy ID
    to either the NCBI taxonomy ID of closest or list of secondary closest evolutionary relatives, or appropriate UCSC
    fasta record handles (genome assembly abbreviations used in record IDs).
    """
    config, taxid_dict = directoryUtility.config, directoryUtility.taxid_dict
    closest_evo_taxids, second_evo_taxids = dict(config['TAXONOMY_CLOSEST_EVOLUTIONARY_RELATIVE']),\
                                            dict(config['TAXONOMY_SECONDARY_EVOLUTIONARY_RELATIVES'])
    closest_evo_taxids = {int(k):int(closest_evo_taxids[k]) for k in closest_evo_taxids}
    second_evo_taxids = {int(k):[int(taxid.strip()) for taxid in second_evo_taxids[k].split(',')] for k in second_evo_taxids}
    closest_evo_handles, second_evo_handles = {},{}
    for taxid in taxid_dict:
        try:
            cid = closest_evo_taxids[taxid]
            sids = second_evo_taxids[taxid]
        except KeyError:
            raise KeyError("Missing evolutionary relative data for taxid {0} in config file".format(taxid))
        c_handle = ucsc_tax_table.loc[ucsc_tax_table['NCBI_taxid']==cid,'fasta_handle'].iloc[0]
        s_handles = [ucsc_tax_table.loc[ucsc_tax_table['NCBI_taxid'] == sid, 'fasta_handle'].iloc[0] for sid in sids]
        closest_evo_handles[taxid],second_evo_handles[taxid] = c_handle,s_handles
    return closest_evo_taxids,closest_evo_handles,second_evo_taxids,second_evo_handles

def evo_relative_alignment(dir_vars,taxid,ucsc_tid,gid,how='all',alt_outpath=""):
    """For a given species taxid, ucsc transcript id, and NCBI gene ID corresponding to orthologs, provides a profile
    alignment of all NCBI ortholog sequences against the subset of UCSC sequences specified with how.
    :param (dict) dir_vars: contains file paths for run
    :param (int) taxid: Taxonomy ID corresponding to species in NCBI analysis set
    :param ucsc_tid: UCSC transcript ID, identifies correct raw UCSC data file.
    :param gid: NCBI human Gene ID used to identify appropriate orthologs file.
    :param how: 'all': use both closest and secondary evo relative data in alignment. 'closest': closest relative only,
    issue warning and use secondary relatives if no record for closest relative. 'secondary': secondary relatives only
    :return:(pandas.Index) evo_records: Record IDs of used close evolutionary relatives for allNCBI alignment
    """
    allNCBI,reorg_parent,ds_id = [dir_vars[label] for label in ["allNCBI_parent","reorg_parent","dataset_identifier"]]
    ucsc_raw_fpath = "{0}/{1}/{2}.fasta".format(reorg_parent,ds_id,ucsc_tid)
    ucsc_raw_df = UCSC_fasta_df(ucsc_raw_fpath)
    orthologs_aa = dir_vars['orthologs_aa']
    ortho_fpath = "{0}/{1}/{2}.fasta".format(orthologs_aa,taxid,gid)
    if not os.path.exists(ortho_fpath):
        emsg = "No NCBI record data for taxid {0} present.".format(taxid)
        raise UCSCerrors.SequenceDataError(0,emsg)
    else:
        #Secondary case for NCBI ortholog Gene ID having no associated Protein IDs (i.e. Gene ID 101721569)
        ortho_ids = [fasta.id for fasta in SeqIO.parse(ortho_fpath, 'fasta')]
        if len(ortho_ids) == 0:
            emsg = "No NCBI record data for taxid {0} present.".format(taxid)
            raise UCSCerrors.SequenceDataError(0, emsg)
    if alt_outpath:
        er_aln_outpath = alt_outpath
    else:
        er_aln_outpath = "{0}/evo_relative_alignments/{1}/{2}.fasta".format(allNCBI, taxid, ucsc_tid)
    #Identify closest evolutionary relatives based on how
    closest_handle, secondary_handles = CLOSEST_EVO_HANDLES[taxid],SECONDARY_EVO_HANDLES[taxid]
    if how == 'all':
        handles = [closest_handle]
        handles.extend(secondary_handles)
    elif how == 'closest':
        handles = [closest_handle]
    elif how == 'secondary':
        handles = secondary_handles
    else:
        raise ValueError("parameter 'how' must be 'all','closest',or 'secondary'.")
    ce_pat = "|".join(handles)
    ce_fpath, all_orthologs_fpath = "tmp/evo_relatives_subset.fasta", "tmp/all_orthologs.fasta"
    evo_records = ucsc_raw_df.index[ucsc_raw_df.index.str.contains(ce_pat)]
    evo_records = evo_records[ucsc_raw_df.loc[evo_records,'length']!=0]
    if len(evo_records) == 0 and how == 'closest':
        #If how==closest but no closest evo data, use secondary evolutionary relative data
        ce_pat = "|".join(secondary_handles)
        ce_fpath, all_orthologs_fpath = "tmp/evo_relatives_subset.fasta", "tmp/all_orthologs.fasta"
        evo_records = ucsc_raw_df.index[ucsc_raw_df.index.str.contains(ce_pat)]
        evo_records = evo_records[ucsc_raw_df.loc[evo_records, 'length'] != 0]
    if len(evo_records) == 0:
        #If all closest evo records are missing, if only one ortholog sequence, write to all_orthologs_fpath. else raise
        #SequenceDataError
        ortho_ids = [fasta.id for fasta in SeqIO.parse(ortho_fpath,'fasta')]
        if len(ortho_ids) == 1:
            filter_fasta_infile(ortho_ids, ortho_fpath, outfile_path=all_orthologs_fpath)
            return []
        else:
            emsg = "No closest evolutionary relative UCSC data, cannot filter multiple NCBI records."
            raise UCSCerrors.SequenceDataError(2,emsg)
    else:
        #, filter closest evo records and aligned orthologs for taxid into tmp files
        filter_fasta_infile(evo_records,ucsc_raw_fpath,outfile_path=ce_fpath)
        fautil.MSA(ortho_fpath,all_orthologs_fpath)
        fautil.profile_MSA(ce_fpath,all_orthologs_fpath,er_aln_outpath)
        return evo_records

def write_bestNCBI_record(dir_vars,taxid,ucsc_tid,NCBI_hgid,evo_record_idx,ucsc_clade='all',how='best_raw'):
    """For a given UCSC transcript ID, corresponding NCBI gene ID, and taxid, selects the best NCBI record by maximum
    identity to records in evo_record_idx. Writes best raw NCBI record file and bestNCBI record against UCSC MSA file.
    :param dir_vars: Contains directory paths. See utility/directoryUtility.py
    :param taxid: Int Taxonomy ID for species
    :param ucsc_tid: UCSC transcript identifier used to name appropriate files
    :param NCBI_gid: NCBI Human Gene ID used to locate ortholog NCBI data
    :param (array-like) evo_record_idx: Contains record ids corresponding to UCSC data for
    :return:
    """
    orthologs_aa,allNCBI,bestNCBI  = [dir_vars[label] for label in ['orthologs_aa','allNCBI_parent','bestNCBI_parent']]
    ucsc_raw_fpath = "{0}/{1}.fasta".format(dir_vars['UCSC_raw_parent'],ucsc_tid)
    ortholog_fasta_fpath = "{0}/{1}/{2}.fasta".format(orthologs_aa,taxid,NCBI_hgid)
    evo_rel_aln_fpath = "{0}/evo_relative_alignments/{1}/{2}.fasta".format(allNCBI,taxid,ucsc_tid)
    ortholog_fasta_ids = [fasta.id for fasta in SeqIO.parse(ortholog_fasta_fpath,'fasta')]
    best_raw_fpath= "{0}/best_raw/{1}/{2}.fasta".format(bestNCBI,taxid,ucsc_tid)
    best_combined_aln_fpath = "{0}/combined/{1}/{2}.fasta".format(bestNCBI,taxid,ucsc_tid)

    if ucsc_clade != 'all':
        clade_generator = fautil.UCSC_subtax_generator(ucsc_raw_fpath,ucsc_clade)
        tmp_ucsc_clade_fpath = "tmp/ucsc_clade.fasta"
        SeqIO.write(clade_generator, tmp_ucsc_clade_fpath, 'fasta')
        ucsc_raw_fpath = tmp_ucsc_clade_fpath

    if len(ortholog_fasta_ids) == 1:
        #Only one ortholog - no min dist record selection, write bestNCBI 1) raw ortholog file and 2) combined alignment
        filter_fasta_infile(ortholog_fasta_ids,ortholog_fasta_fpath,best_raw_fpath)
        bestNCBI_record_id = ortholog_fasta_ids[0]
    else:
        evo_aln_records_df = fautil.load_UCSC_NCBI_df(evo_rel_aln_fpath)
        evo_aln_fasta_ids = evo_aln_records_df.index
        id_dm, align_srs = construct_id_dm(evo_aln_fasta_ids, evo_rel_aln_fpath,filtered=True,aligned=True)
        md_row, md = min_dist_spec_record(id_dm,align_srs.index,ortholog_fasta_ids,evo_record_idx,evo_aln_records_df)
        bestNCBI_record_id= md_row.name
        filter_fasta_infile([bestNCBI_record_id], ortholog_fasta_fpath, best_raw_fpath)

    if how == 'best_raw':
        fautil.profile_MSA(ucsc_raw_fpath, best_raw_fpath, best_combined_aln_fpath)
    elif how == 'best_evo_realn':
        tmp_evo_records, tmp_bestevo_realn = "tmp/evo_records.fasta", "tmp/best_evo_aln_realigned.fasta"
        filter_fasta_infile(evo_record_idx, ucsc_raw_fpath, outfile_path=tmp_evo_records)
        fautil.profile_MSA(tmp_evo_records, best_raw_fpath, out_fpath=tmp_bestevo_realn)
        fautil.profile_MSA(ucsc_raw_fpath, tmp_bestevo_realn, best_combined_aln_fpath)
    return bestNCBI_record_id

def length_metrics_df(lm_fpath, taxid_dict):
    """Reads or initializes a DataFrame from lm_fpath. Columns are based on evolutionary relative data read from config
    file.
    :param lm_fpath: File path to length_metrics.tsv. If no file exists at path, makes new DataFrame
    :param taxid_dict: Contains Taxonomy IDs of all NCBI species for analysis as keys.
    :return: lm_df: DataFrame indexed on UCSC transcript IDs, contains NCBI Human Gene ID and lengths of sequences for
    UCSC evolutionary relatives and NCBI records.
    """
    core_labels = ["clade_mean", "clade_median", "bestncbi_mean", "bestncbi_median"]
    closest_evo_taxids = list(CLOSEST_EVO_TAXIDS.values())
    second_evo_taxids = []
    for sids in SECONDARY_EVO_TAXIDS.values():
        second_evo_taxids.extend(sids)
    ucsc_taxids = closest_evo_taxids + second_evo_taxids
    ordered = OrderedDict(zip(ucsc_taxids, ucsc_taxids))
    ordered_ucsc_taxids = list(ordered.keys())
    ucsc_labels = ["{0}_length".format(taxid) for taxid in ordered_ucsc_taxids]
    ncbi_labels = ["{0}_ncbi_length".format(taxid) for taxid in taxid_dict]
    lm_columns = ['NCBI_hgid']+core_labels + ucsc_labels + ncbi_labels
    if not os.path.exists(lm_fpath):
        lm_df = pd.DataFrame(columns=lm_columns,dtype=float)
        lm_df['NCBI_hgid'] = lm_df['NCBI_hgid'].astype('str')
        lm_df.index.name = "UCSC_transcript_id"
    else:
        field_types = ['str']
        field_types.extend(['float'] * (len(lm_columns)-1))
        field_conv_dict = dict(zip(lm_columns, field_types))
        lm_df = pd.read_csv(lm_fpath, sep='\t', index_col="UCSC_transcript_id", dtype=field_conv_dict)
    return lm_df, ordered_ucsc_taxids

def length_metrics_row(dir_vars,taxid_dict,lm_df,ordered_ucsc_taxids,ucsc_tid,ncbi_hgid,ucsc_clade='all'):
    """Updates lm_df with length information on clade, evolutionary relatives, and ncbi records. Should only be used
    after best record selection/ write_bestNCBI_record has been called for each species in taxid_dict.
    :param (dict) dir_vars: see utility/directoryUtility.py
    :param (dict) taxid_dict: see utility/directoryUtility.py
    :param lm_df: DataFrame with columns as defined in length_metrics_df.
    :param ordered_ucsc_taxids: Ordered, non-redundant collection of UCSC taxids of evolutionary relatives for species
    in NCBI analysis set; returned by length_metrics_df
    :param ucsc_tid: Transcript ID used to locate raw UCSC files
    :param ncbi_hgid: NCBI Human Gene ID correspondgin to ucsc_tid
    :param ucsc_clade: identifier corresponding to subgroup of UCSC records. See utility/fastaUtility.py and
    UCSC_clade_positions for details.
    :return: lm_df, updated with length information for this transcript id and gene id.
    """
    reorg_parent, ds_id = [dir_vars[label] for label in ['reorg_parent','dataset_identifier']]
    ucsc_raw_fpath = "{0}/{1}/{2}.fasta".format(reorg_parent,ds_id,ucsc_tid)
    ucsc_df = fautil.UCSC_fasta_df(ucsc_raw_fpath)
    for taxid in ordered_ucsc_taxids:
        ucsc_col_label = "{0}_length".format(taxid)
        tax_len = ucsc_df.loc[(ucsc_df['NCBI_taxid'] == taxid), 'length'].iloc[0]
        lm_df.loc[ucsc_tid,ucsc_col_label] = tax_len
    clade_df, rest = fautil.partition_UCSC_by_clade(ucsc_df,clade=ucsc_clade)
    lm_df.loc[ucsc_tid, 'clade_mean'] = np.mean(clade_df['length'])
    lm_df.loc[ucsc_tid, 'clade_median'] = np.median(clade_df['length'])
    lm_df.loc[ucsc_tid,'NCBI_hgid'] = str(ncbi_hgid)
    ncbi_lens = []
    for taxid in taxid_dict:
        bestNCBI_fpath = "{0}/best_raw/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'],taxid,ucsc_tid)
        if os.path.exists(bestNCBI_fpath):
            ncbi_df = fautil.NCBI_fasta_df(bestNCBI_fpath,taxid_dict)
            # tax_len = ncbi_df.loc[(ncbi_df['NCBI_taxid'] == taxid), 'length'].iloc[0]
            tax_len = ncbi_df["length"].iloc[0]
            ncbi_col_label = "{0}_ncbi_length".format(taxid)
            lm_df.loc[ucsc_tid, ncbi_col_label] = tax_len
            ncbi_lens.append(tax_len)
    lm_df.loc[ucsc_tid,'bestncbi_mean'] = np.mean(ncbi_lens)
    lm_df.loc[ucsc_tid, 'bestncbi_median'] = np.median(ncbi_lens)
    return lm_df

def filter_ortholog_data(force_overwrite_tid_subset=[],clade="all",alt_lm_fpath=""):
    """Master function for carrying out allNCBI/ bestNCBI file writing. In allNCBI: writes all ortholog against closest
    evolutionary relative alignments (or skips if no evolutionary relative data and multiple NCBI ortholog sequences).
    In bestNCBI: writes 1) unaligned best identity record to evolutionary relative data and 2) combined ucsc
    clade-filtered against best NCBI record profile alignment.
    :param force_overwrite_tid_subset: array-like, if provided forces alignment file rewriting and record length
    recalculations for UCSC transcript IDs in it
    :param clade: Must correspond to accepted values in fastaUtility ('all','mammalia', 'boreoeutheria',
    'euarchontoglires', 'primates').
    :param alt_lm_fpath: If provided, reads and writes length_metrics table to alternative file path.
    :return:
    """
    from utility.directoryUtility import taxid_dict,dir_vars
    xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars['xref_summary'])
    xref_table = orutil.load_NCBI_xref_table(xref_fpath)
    NCBI_gid_table = xref_table.loc[:, ["stable_gene_id", "NCBI_gid", "HGNC_symbol"]]
    if alt_lm_fpath:
        lm_fpath = alt_lm_fpath
    else:
        lm_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
    lm_df, ordered_ucsc_taxids = length_metrics_df(lm_fpath,taxid_dict)

    orthologs_fpath = "{0}/orthologs_final.tsv".format(dir_vars['orthologs_parent'])
    ortholog_id_table = orutil.load_orthologs_table(orthologs_fpath)
    non_empty_orthoIDs = ortholog_id_table.dropna(how='all')
    some_NCBI_IDs = NCBI_gid_table.loc[(NCBI_gid_table["NCBI_gid"].isin(non_empty_orthoIDs.index.astype(str))), :]
    fpaths, check_flags, error_dfs = UCSCerrors.load_all_taxid_errors(dir_vars, taxid_dict, log_type='record',
                                                                      error_type="SequenceDataError",error_codes=[0,2])
    for ucsc_tid, row in some_NCBI_IDs.iterrows():
        ncbi_hgid = row["NCBI_gid"]
        for taxid in taxid_dict:
            combined_best_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'],taxid,ucsc_tid)
            error_fpath,error_check,error_df = fpaths[taxid], check_flags[taxid], error_dfs[taxid]
            if (not os.path.exists(combined_best_fpath)) or ucsc_tid in force_overwrite_tid_subset \
                    or (ucsc_tid in lm_df.index and lm_df.loc[ucsc_tid, :].isnull().all()):
                if error_check and ucsc_tid in error_df['tid'].unique():
                    UCSCerrors.print_errors(error_df,ucsc_tid)
                    continue
                try:
                    evo_records = evo_relative_alignment(dir_vars,taxid,ucsc_tid,ncbi_hgid)
                    write_bestNCBI_record(dir_vars,taxid,ucsc_tid,ncbi_hgid,evo_records,ucsc_clade=clade,how='best_raw')
                except UCSCerrors.SequenceDataError as sde:
                    UCSCerrors.write_errors(error_fpath,ucsc_tid,sde)
                    continue
        if ucsc_tid not in lm_df.index or ucsc_tid in force_overwrite_tid_subset :
            lm_df = length_metrics_row(dir_vars,taxid_dict,lm_df,ordered_ucsc_taxids,ucsc_tid,ncbi_hgid,ucsc_clade=clade)
            lm_df.to_csv(lm_fpath, sep='\t', float_format='%.4f')
    lm_df.to_csv(lm_fpath, sep='\t', float_format='%.4f')

#Read evolutionary relative data from config file.
CLOSEST_EVO_TAXIDS,CLOSEST_EVO_HANDLES,SECONDARY_EVO_TAXIDS,SECONDARY_EVO_HANDLES = read_evo_relative_data()

def main():
    pd.options.display.max_columns = None
    filter_ortholog_data(clade="boreoeutheria")

if __name__ == '__main__':
    main()




