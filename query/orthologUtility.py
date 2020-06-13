import pandas as pd
import numpy as np
import os

def load_orthologs_table(orthologs_fpath,gid_dtype='str',index_type='int'):
    ortholog_table = pd.read_csv(orthologs_fpath,sep='\t',index_col=0,dtype=gid_dtype)
    ortholog_table.index = ortholog_table.index.astype(index_type)
    return ortholog_table

def load_NCBI_xref_table(NCBI_xrefs_fpath,gid_dtype='str'):
    xrefs = pd.read_csv(NCBI_xrefs_fpath, dtype=str, delimiter='\t', index_col='UCSC_transcript_id')
    xrefs.loc[:,'NCBI_gid'] = xrefs["NCBI_gid"].astype(gid_dtype)
    return xrefs

def boolean_df_agg(ortholog_df,bool_df,how='any',comp=False):
    if how=='all':
        match_idx = bool_df.all(axis=1)
    elif how=='any':
        match_idx = bool_df.any(axis=1)
    else:
        raise ValueError("value of how must be 'any' or 'all")
    if comp:
        match_df = ortholog_df.loc[~match_idx,:].sort_index()
        return match_df
    else:
        match_df = ortholog_df.loc[match_idx,:].sort_index()
        return match_df

def ortholog_patmatch(ortholog_table,match_pat,how="any",comp=False):
    """Returns DataFrame of rows from orthologs

    :param ortholog_fpath: path to orthologs tsv table. Rows with all null values will be dropped.
    :param pat: re pattern to match against
    :param how: "any" or "all", uses appropriate boolean Series filter method
    :param comp: boolean. If True, returns complement of match_df.index from ortholog_fpath (all null value rows
    will be absent from complement).
    how=any,comp=False: rows with any matches
    how=all,comp=False: rows with all matches
    how=any,comp=True: rows with all mismatches
    how=all,comp=True will return rows with any mismatches
    :return: match_df: DataFrame containing rows from ortholog_fpath for which any entry matches pat.
    """
    from IPython.display import display
    ortholog_table = ortholog_table.dropna(how='all')
    print("Ortholog entries: {0}".format(len(ortholog_table)))
    #Generate boolean series
    bool_df = pd.DataFrame(index=ortholog_table.index, columns=ortholog_table.columns)
    for col in ortholog_table:
        ortholog_col = ortholog_table[col]
        bool_col = ortholog_col.str.contains(match_pat)
        bool_df.loc[:,col] = bool_col
    match_df = boolean_df_agg(ortholog_table,bool_df,how=how,comp=comp)
    return match_df

def ortholog_na_filter(ortholog_table,how="any",comp=False):
    """

    :param ortholog_fpath: path to orthologs tsv table. Rows with all null values will be dropped.
    :param how: "any" or "all", uses appropriate boolean Series filter method
    :param comp: boolean. If True, returns complement of match_df.index from ortholog_fpath
    how=any,comp=False: rows with any nan values
    how=all,comp=False: rows with all nan
    how=any,comp=True: rows with all non-nan
    how=all,comp=True will return rows with any non-nan
    :return: match_df: DataFrame containing rows from ortholog_fpath for which any entry matches pat.
    """
    print("Ortholog entries: {0}".format(len(ortholog_table)))
    bool_df = pd.DataFrame(index=ortholog_table.index, columns=ortholog_table.columns)
    for col in ortholog_table:
        ortholog_col = ortholog_table[col]
        bool_col = ortholog_col.isna()
        bool_df.loc[:, col] = bool_col
    match_df = boolean_df_agg(ortholog_table, bool_df, how=how, comp=comp)
    return match_df

def any_NCBI_xref(alt_xref_fpath="",alt_ortho_fpath=""):
    """Returns DataFrame indexed on UCSC transcript ID containing Ensembl stable GeneID, NCBI GID, and HGNC symbol
    information for all transcript IDs for which at least one NCBI ortholog

    :param alt_xref_fpath:
    :param alt_ortho_fpath:
    :return:
    """
    from utility.directoryUtility import dir_vars
    if alt_xref_fpath:
        xref_fpath = alt_xref_fpath
    else:
        xref_fpath = "{0}/NCBI_xrefs.tsv".format(dir_vars['xref_summary'])
    if alt_ortho_fpath:
        orthologs_fpath = alt_ortho_fpath
    else:
        orthologs_fpath = "{0}/orthologs_final.tsv".format(dir_vars['orthologs_parent'])
    xref_table = load_NCBI_xref_table(xref_fpath)
    NCBI_gid_table = xref_table.loc[:, ["stable_gene_id", "NCBI_gid", "HGNC_symbol"]]

    ortholog_id_table = load_orthologs_table(orthologs_fpath)
    non_empty_orthoIDs = ortholog_id_table.dropna(how='all')
    any_NCBI_IDs = NCBI_gid_table.loc[(NCBI_gid_table["NCBI_gid"].isin(non_empty_orthoIDs.index.astype(str))), :]
    return any_NCBI_IDs