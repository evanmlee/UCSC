import os
import pandas as pd
import numpy as np
import re

from query import sql_orthologs as so,api_orthologs as ao
from query import orthologUtility as orutil
from utility.UCSCerrors import load_errors, NCBIQueryError


def compile_ortholog_data(local_orthologs_fpath, API_orthologs_fpath,tax_id_columns):
    """Combines local and API ortholog data into final table from NCBI gid to analysis species ortholog Gene IDs.
    
    :param local_orthologs_fpath: File path to ortholog table for analysis derived from SQL (NCBI orthologs FTP) data
    :param API_orthologs_fpath: File path to ortholog table for analysis derived from queries to NCBI/Gene website 
    (i.e. for gene ids which were absent in local_orthologs or for species absent from SQL data) 
    :param (dict) tax_id_columns: keys - NCBI Tax IDs for species for which NCBI seq data is being used; values: column
    labels in ortholog tables (generally {taxid}_gid)
    :return: final orthologs table
    Final orthologs table is filled in preferential order: 1) local orthologs data for species for which NCBI orthologs
    FTP file is available 2) API orthologs data for remaining Gene IDs not represented in FTP data or for all Gene IDs
    if species is not present in FTP data.
    contains
    """
    local_orthologs = pd.read_csv(local_orthologs_fpath,sep='\t',dtype='str',index_col=0)
    local_entry_counts = local_orthologs.count()
    local_missing_cols = [col for col,count in local_entry_counts.iteritems() if count==0]
    API_orthologs = pd.read_csv(API_orthologs_fpath, sep='\t', dtype='str', index_col=0)
    union_index = local_orthologs.index.union(API_orthologs.index)
    final_orthologs = pd.DataFrame(index=union_index,columns=API_orthologs.columns)
    final_orthologs.loc[local_orthologs.index,:] = local_orthologs
    API_local_comp = API_orthologs.index.difference(local_orthologs.index)
    local_API_int = local_orthologs.index.intersection(API_orthologs.index)
    final_orthologs.loc[API_local_comp,:] = API_orthologs.loc[API_local_comp,:]
    final_orthologs.loc[local_API_int,local_missing_cols] = API_orthologs.loc[local_API_int,local_missing_cols]
    return final_orthologs

def _compare_formatted_orthologs(suspect_fpath,correct_fpath):
    """Internal function. Compares entries from two versions of the formatted_ortholgos table returned by
    format_orthologs_sql. Returns a DataFrame indexed on NCBI human gid and containing rows for which entries
    differ between the tables at suspect_fpath and correct_fpath.

    :param suspect_fpath: file path to formatted orthologs table with suspected incorrect entries
    :param correct_fpath: file path to formatted_orthologs table with correct entries (ie populated locally using
    gene_orthologs SQL table).
    NOTE: columns must be in same order in both dataframes.
    :return: diff_df: DataFrame object, indexed on NCBI hgid, contains rows from correct_fpath table for which any
    entry is different in suspect_fpath (ie missing values, incorrect GIDs).
    """
    suspect_df = pd.read_csv(suspect_fpath,sep='\t',dtype='str',index_col=0)
    correct_df = pd.read_csv(correct_fpath, sep='\t', dtype='str',index_col=0)
    suspect_df.index.name = correct_df.index.name
    suspect_df.columns = correct_df.columns
    diff_df = pd.DataFrame(columns=correct_df.columns)
    for hgid, correct_row in correct_df.iterrows():
        if hgid not in suspect_df.index:
            diff_df.loc[hgid, :] = correct_row
            continue
        suspect_row = suspect_df.loc[hgid,:]
        for col,entry in correct_row.iteritems():
            sus_entry = suspect_row[col]
            if (type(sus_entry) == str and type(entry) == str and
                sus_entry != entry) or \
                (type(sus_entry) == float and np.isnan(sus_entry) and type(entry) == str) or \
                (type(entry) == float and np.isnan(entry) and type(sus_entry) == str):
                diff_df.loc[hgid,:] = correct_row
    return diff_df

def _identify_suspect_rows(local_orthologs_fpath,API_orthologs_fpath,match_pat,include_na=True):
    """

    :param local_orthologs_fpath:
    :param API_orthologs_fpath:
    :param include_na:
    :return: Rows from API orthologs not present in local orthologs which have pattern mismatches to match_pat
    """
    match_pat = "^10|^11"
    API_orthologs = orutil.load_orthologs_table(API_orthologs_fpath)
    print("API orthologs: {0}".format(len(API_orthologs)))
    local_orthologs = orutil.load_orthologs_table(local_orthologs_fpath)
    print("local orthologs: {0}".format(len(local_orthologs)))
    any_mismatch = orutil.ortholog_patmatch(API_orthologs,match_pat,how='all',comp=True)
    # API_all_valid = API_orthologs.dropna(how='any')
    if include_na:
        any_na = orutil.ortholog_na_filter(API_orthologs,how='any',comp=False)
        suspect_idx = any_na.index.union(any_mismatch.index)
        nla_suspect_fpath = "tmp/nonlocal_API_suspect.tsv"
    else:
        suspect_idx = any_mismatch.index
        nla_suspect_fpath = "tmp/nonlocal_API_suspect_mismatch_only.tsv"
    #Identify records in API orthologs (longer list) not in locally (SQL) matched orthologs
    non_local_API_idx = API_orthologs.index.difference(local_orthologs.index)
    non_local_API = API_orthologs.loc[non_local_API_idx,:]
    non_local_API_fpath = "tmp/nonlocal_API.tsv"
    non_local_API.to_csv(non_local_API_fpath,sep='\t')

    nla_suspect_idx = non_local_API.index.intersection(suspect_idx)
    non_local_API_suspect = non_local_API.loc[nla_suspect_idx,:]
    non_local_API_suspect.to_csv(nla_suspect_fpath,sep='\t')
    return non_local_API_suspect

def missing_local_orthologs(orthologs_table,taxid_cols,local_orthologs):
    """Checkes orthologs_table for rows missing value for specified taxid, only returning rows present in
    local_orthologs.indez

    :param (DataFrame) orthologs_table: DataFrame indexed on NCBI human gid and containing columns named [taxid]_gid.
     Required column correspondign to taxid
    :param (dict) taxid_cols: Maps taxids to column names for species which are missing from local_orthologs
    :param (DataFrame) local_orthologs: DataFrame corresponding to mapped orthologs pulled from SQL gene_orthologs table
    :return: DataFrame containing subset of orthologs table with missing values in columns speciefied by taxid_cols
    and whose index is in local_orthologs.
    """
    sub_table = orthologs_table.loc[local_orthologs.index,taxid_cols.values()]
    bool_df = sub_table.isna()
    return orutil.boolean_df_agg(sub_table,bool_df,how='any',comp=False)


def write_ortholog_errors(taxid_dict,dir_vars,outpath="",overwrite=False):
    from IPython.display import display
    taxid_str = ",".join(taxid_dict.keys())

    xref_summary,orthologs_dir, summary_run = [dir_vars[label] for label in  ['xref_summary','orthologs_parent',
                                                                              'summary_run']]
    local_orthologs_fpath = "{0}/knownCanonical_orthologs_local.tsv".format(orthologs_dir)
    xref_fpath = "{0}/NCBI_xrefs.tsv".format(xref_summary)
    # final_errors_fpath = "{0}/errors.tsv".format(summary_run)
    final_orthologs_fpath = "{0}/orthologs_final.tsv".format(orthologs_dir)
    query_errors_fpath = "{0}/errors_ncbiquery.tsv".format(summary_run)
    length_metrics_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
    if not outpath:
        outpath = "{0}/ortholog_errors.tsv".format(summary_run)

    xref_df = pd.read_csv(xref_fpath,sep='\t',index_col=0)
    local_orthologs_df = pd.read_csv(local_orthologs_fpath, sep='\t', index_col=0)

    local_allna = orutil.ortholog_na_filter(local_orthologs_df,how='all')
    final_orthologs_df = pd.read_csv(final_orthologs_fpath,sep='\t',index_col=0)
    final_allna = orutil.ortholog_na_filter(final_orthologs_df,how='all')
    lm_df = pd.read_csv(length_metrics_fpath,sep='\t',index_col=0)

    check_qes, query_errors_df = load_errors(query_errors_fpath)
    qe_gids = xref_df.loc[query_errors_df['tid'],'NCBI_gid']

    error_columns = ['tid', 'error_type', 'error_code', 'error_message']
    ortholog_errors = pd.DataFrame(columns=error_columns)
    if not os.path.exists(outpath) or overwrite:
        for transcript_id, row in xref_df.iterrows():
            hgid = row['NCBI_gid']
            if transcript_id in local_allna.index:
                ncbi_qe = NCBIQueryError(1,"No SQL orthologs IDs available for hgid {0}".format(ncbi_qe))
            elif transcript_id in query_errors_df['tid'].unique():
                ncbi_qe = NCBIQueryError(0, "No orthologs available for NCBI hgid {0}".format(hgid))
            elif hgid in qe_gids.unique():
                ncbi_qe = NCBIQueryError(0, "No orthologs available for NCBI hgid {0}".format(hgid))
            elif hgid in final_allna.index:
                ncbi_qe = NCBIQueryError(2, "Tax IDs {0} had no orthologs available for NCBI hgid {1}".format(taxid_str,hgid))
            elif transcript_id not in lm_df.index:
                ncbi_qe = NCBIQueryError(3, "Transcript ID dropped from final analysis subset")
            else:
                continue
            etype, ecode, emsg = ncbi_qe.error_type, ncbi_qe.code, ncbi_qe.message
            error_row = pd.Series({'tid':transcript_id,'error_type': etype, 'error_code': ecode, 'error_message': emsg})
            ortholog_errors = ortholog_errors.append(error_row,ignore_index=True)

        ortholog_errors.to_csv(outpath,sep='\t')
    else:
        ortholog_errors = pd.read_csv(outpath,sep='\t',index_col=0)

    check_paths = True
    if check_paths:
        ambig_errors = ortholog_errors.loc[ortholog_errors['error_code']==3,:]
        ambig_gids = xref_df.loc[ambig_errors['tid'],'NCBI_gid']
        dd_gids = ambig_gids.drop_duplicates()
        ambig_orthologs = final_orthologs_df.loc[dd_gids,:]
        display(ambig_orthologs)

        orthologs_aa, allncbi_dir = dir_vars['orthologs_aa'],dir_vars['allNCBI_parent']
        ambig_orthologs_df_count = ambig_orthologs.count()
        print("Valid ortholog gids for ambiguous error GIDs: {0}".format(sum(ambig_orthologs_df_count)))

        file_exists_count = 0
        for gid in ambig_orthologs.index:
            for taxid in taxid_dict:
                ortho_fpath = "{0}/{1}/{2}.fasta".format(orthologs_aa,taxid,gid)
                if os.path.exists(ortho_fpath):
                    file_exists_count += 1
        print("Existing ortholog file path count for ambigous error orthologs: {0}".format(file_exists_count))

        allncbi_fpath_count = 0
        for tid in ambig_errors['tid']:
            allncbi_fpath = "{0}/combined/{1}.fasta".format(allncbi_dir,tid)
            if os.path.exists(allncbi_fpath):
                allncbi_fpath_count+=1
        print("Number of UCSC Transcript IDs with 'dropped from analysis' errors: {0}".format(len(ambig_errors)))
        print("Existing allNCBI file path count for ambigous error tids: {0}".format(allncbi_fpath_count))

def allseq_NCBI_UCSC_slignment(NCBI_xref_df, taxid_dict, dir_vars,gid_subset=[]):
    """File preparation for allNCBI data. Prepares: 1) compiled NCBI raw and aligned data and 2) all NCBI data profiled
    aligned against corresponding UCSC data.

    :param NCBI_xref_df:
    :param taxid_dict:
    :param dir_vars:
    :param gid_subset:
    :return:
    """
    from utility.fastaUtility import MSA_fpath_list, profile_MSA
    # Directory paths
    dir_labels = ["orthologs_parent","orthologs_aa","UCSC_raw_parent","allNCBI_parent"]
    orthologs_dir, orthologs_aa, ucsc_raw, allNCBI_parent = [dir_vars[label] for label in dir_labels]
    subdirs = ["NCBI_raw","NCBI_alignments","combined"]
    allNCBI_raw,allNCBI_alignments,allNCBI_combined = ["{0}/{1}".format(allNCBI_parent,sub) for sub in subdirs]

    for UCSC_tid, row in NCBI_xref_df.iterrows():
        NCBI_gid = row["NCBI_gid"]
        force_filewrite = False
        if len(gid_subset) > 0:
            if NCBI_gid not in gid_subset:
                continue
            else:
                force_filewrite = True
        fpath_list = []
        for taxid in taxid_dict:
            NCBI_fpath = "{0}/{1}/{2}.fasta".format(orthologs_aa, taxid, NCBI_gid)
            if (os.path.exists(NCBI_fpath)):
                fpath_list.append(NCBI_fpath)
        if len(fpath_list) > 0:
            NCBI_raw_fpath = "{0}/{1}.fasta".format(allNCBI_raw, NCBI_gid)
            NCBI_alignment_fpath = "{0}/{1}.fasta".format(allNCBI_alignments, NCBI_gid)
            if not os.path.exists(NCBI_alignment_fpath) or force_filewrite:
                MSA_fpath_list(fpath_list, NCBI_raw_fpath, NCBI_alignment_fpath)

            combined_alignment_fpath = "{0}/{1}.fasta".format(allNCBI_combined, UCSC_tid)
            if not os.path.exists(combined_alignment_fpath) or force_filewrite:
                UCSC_raw_tid_fasta = "{0}/{1}.fasta".format(ucsc_raw, UCSC_tid)
                try:
                    assert (os.path.exists(UCSC_raw_tid_fasta))
                    assert (os.path.exists(NCBI_alignment_fpath))
                    profile_MSA(in1_fpath=UCSC_raw_tid_fasta, in2_fpath=NCBI_alignment_fpath,
                                out_fpath=combined_alignment_fpath)
                except AssertionError as e:
                    print("UCSC_tid_fasta_path: {0}".format(UCSC_raw_tid_fasta))
                    print("NCBI_alignment_fpath: {0}".format(NCBI_alignment_fpath))

def main():
    from IPython.display import display
    from utility.directoryUtility import config,taxid_dict,dir_vars
    orthologs_dir = dir_vars["orthologs_parent"]
    reorg_dir, dataset_name,summary_run_dir = dir_vars['reorg_parent'],dir_vars['dataset_name'],dir_vars['summary_run']

    formatted_local_fpath = "{0}/knownCanonical_orthologs_local.tsv".format(orthologs_dir)
    local_orthologs_fpath = formatted_local_fpath
    API_orthologs_fpath = "{0}/NCBI_orthologs_API.tsv".format(orthologs_dir)
    NCBI_xref_inpath = "{0}/xref_summary_{1}/NCBI_xrefs.tsv".format(reorg_dir,dataset_name)
    dataset_id = dir_vars["dataset_identifier"]
    final_errors_fpath = "{0}/errors.tsv".format(summary_run_dir)
    final_orthologs_fpath = "{0}/orthologs_final.tsv".format(orthologs_dir)
    sus_fpath = "{0}/NCBI_orthologs_031320.tsv".format(orthologs_dir)

    format_sql_check = False
    if format_sql_check:
        # format orthologs table
        formatted_df, absent_taxids = so.format_orthologs_sql(dir_vars, taxid_dict, formatted_local_fpath)
        display(formatted_df)

    id_mismatch_gids = False
    if id_mismatch_gids:
        valid_match_pat = "^10|^11"
        sus_ortho_df = orutil.load_orthologs_table()
        mismatch_df = orutil.ortholog_patmatch(sus_ortho_df, valid_match_pat, how='all', comp=True)
        with pd.option_context('display.max_columns', None):
            display(mismatch_df)
        suspect_entries_fpath = "{0}/NCBI_orthologs_suspect.tsv".format(orthologs_dir)
        mismatch_df.to_csv(suspect_entries_fpath, sep='\t')

    compile_final_orthologs = True
    if compile_final_orthologs:
        taxid_columns = dict(zip(taxid_dict.keys(),['{0}_gid'.format(k) for k in taxid_dict.keys()]))
        final_orthologs = compile_ortholog_data(local_orthologs_fpath,API_orthologs_fpath,taxid_columns)
        # final_orthologs.to_csv(final_orthologs_fpath,sep='\t')

    check_diff = False
    if check_diff:
        # diff_df = _compare_formatted_orthologs(sus_fpath,API_orthologs_fpath)
        diff_df = _compare_formatted_orthologs(sus_fpath,final_orthologs_fpath)
        display(diff_df)
        diff_df.to_csv('tmp/diff.tsv',sep='\t')

    # Optional repeat geneID mapping to fix na entries or just mismatch entries (ie invalid GeneID values).
    repeat_na_ID_mapping = False
    if repeat_na_ID_mapping:
        non_local_API_suspect = _identify_suspect_rows(local_orthologs_fpath,API_orthologs_fpath,valid_match_pat,include_na=True)
        ao.map_AGS_geneIDs(NCBI_xref_inpath,taxid_dict,API_orthologs_fpath,final_errors_fpath,overwrite_gid=non_local_API_suspect.index)

    repeat_mismatch_mapping = False
    if repeat_mismatch_mapping:
        final_any_mismatch = orutil.ortholog_patmatch(final_orthologs,valid_match_pat,how='all',comp=True)
        ao.map_AGS_geneIDs(NCBI_xref_inpath,taxid_dict,API_orthologs_fpath,final_errors_fpath,overwrite_gid=final_any_mismatch.index)

    # Update protein records for diff_df
    update_pids = False
    if update_pids:
        api_key = config["DIRECTORY"]["NCBIAPIKey"]
        pid_df = ao.fetch_NCBI_protein_records(final_orthologs,taxid_dict,dir_vars,
                                        overwrite_fasta=diff_df.index,NCBI_api_key=api_key)
    update_errors = True
    if update_errors:
        query_error_path = "{0}/errors_ncbiquery.tsv".format(dir_vars['summary_run'])
        outpath = "{0}/ortholog_errors.tsv".format(dir_vars['summary_run'])
        write_ortholog_errors(taxid_dict,dir_vars,outpath=outpath)

    process_allNCBI = False
    if process_allNCBI:
        NCBI_xrefs = orutil.load_NCBI_xref_table(NCBI_xref_inpath, gid_dtype='int')
        allseq_NCBI_UCSC_slignment(NCBI_xrefs,taxid_dict,dir_vars,gid_subset=diff_df.index.values)

if __name__ == "__main__":
    main()

