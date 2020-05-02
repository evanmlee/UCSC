#TODO: Clean up, debug
import os
import pandas as pd
import numpy as np
import re
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException,NoSuchElementException
from UCSCerrors import DataFileError, NCBIQueryError
from UCSCerrors  import write_errors, load_errors
import subprocess
import requests
import xml.etree.ElementTree as ET

def headless_driver():
    #Headless driver with default options
    from selenium.webdriver.chrome.options import Options
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    WINDOW_SIZE = "1920,1080"
    chrome_options.add_argument("--window-size=%s" % WINDOW_SIZE)
    driver = webdriver.Chrome(options=chrome_options)
    return driver

def load_orthologs_table(orthologs_fpath):
    ortholog_table = pd.read_csv(orthologs_fpath,sep='\t',index_col=0,dtype=str)
    # ortholog_table.index = ortholog_table.index.astype('str')
    return ortholog_table

def load_NCBI_xref_table(NCBI_xrefs_fpath,gid_dtype='str'):
    xrefs = pd.read_csv(NCBI_xrefs_fpath, dtype=str, delimiter='\t', index_col='UCSC_transcript_id')
    xrefs.loc[:,'NCBI_gid'] = xrefs["NCBI_gid"].astype(gid_dtype)
    return xrefs

def SQL_transcript_table(dir_vars,tid_index,write_table=False):
    """Returns and optionally writes to file an int-indexed table of UCSC transcript IDs, and corresponding Ensembl
    transcript stable IDs and versions.

    :param (dict) dir_vars: contains path variables for run directories
    :param (pandas Index) tid_index: Index of UCSC_transcript_ids (ENST############.##)
    :param (boolean) write_table: If true, write to sql_transcript_fpath
    :return: DataFrame containing int-indexed UCSC_transcript IDs as well as stable transcript ID and version numbers
    which will be used to query Ensembl SQL tables.
    """
    reorg_parent,dataset_name = dir_vars["reorg_parent"],dir_vars["dataset_name"]
    sql_transcript_fpath = "{0}/{1}_transcript_table.tsv".format(reorg_parent, dataset_name)
    sql_transcript_table_index = tid_index.to_series()
    sql_transcript_table = pd.DataFrame(index=sql_transcript_table_index,
                                        columns=['stable_transcript_id', 'version'], dtype='str')
    sql_transcript_table.loc[:, 'stable_transcript_id'] = sql_transcript_table_index.str.extract("(ENST\d+)\.\d+")
    sql_transcript_table.loc[:, 'version'] = sql_transcript_table_index.str.extract("ENST\d+\.(\d+)")
    int_indexed = sql_transcript_table.reset_index(drop=False)
    int_indexed.index.name = "UCSC_transcript_index"
    if write_table:
        int_indexed.to_csv(sql_transcript_fpath,sep='\t',header=False)
    return int_indexed

def format_orthologs_sql(dir_vars,taxid_dict,formatted_outpath):
    """Reformats NCBI ortholog gid information from the SQL database tables to store appropriate orthologs in separate
    columns for each species.

    :param dir_vars: dictionary of directory paths based on run parameters, see directoryUtility
    :param taxid_dict: dictionary mapping NCBI taxids of species for which NCBI data will be fetched to species name
    :return: formatted_df: DataFrame indexed on NCBI human gid, columns named to [taxid]_gid and populated with
            NCBI gids of orthologs from appropriate species, based on data present in orthologs_sql_fpath.
    :return: absent_taxids: list of taxids that are in taxid_dict but not present in the ortholog_sql table
    """
    from IPython.display import display
    orthologs_parent,dataset_name = dir_vars["orthologs_parent"],dir_vars["dataset_name"]
    orthologs_sql_fpath = "{0}/{1}_orthologs_sql.tsv".format(orthologs_parent,dataset_name)
    if not os.path.exists(orthologs_sql_fpath):
        raise DataFileError(0,"Missing run orthologs file at path: {0}".format(orthologs_sql_fpath))
    col_names = ['ortholog_id','tax_id','GeneID','relationship','Other_tax_id','Other_GeneID']
    orthologs_sql_df = pd.read_csv(orthologs_sql_fpath,sep='\t',dtype=str,index_col=0,names=col_names)
    NCBI_gid_index = pd.Index(orthologs_sql_df["GeneID"].drop_duplicates())
    formatted_columns = ["{0}_gid".format(taxid) for taxid in taxid_dict]
    formatted_df = pd.DataFrame(index=NCBI_gid_index,columns=formatted_columns)
    for _,row in orthologs_sql_df.iterrows():
        hgid = row["GeneID"]
        spec_tax_id = row["Other_tax_id"]
        spec_gid = row["Other_GeneID"]
        spec_col = "{0}_gid".format(spec_tax_id)
        formatted_df.loc[hgid,spec_col] = spec_gid
    formatted_df.index.name = "NCBI_hgid"
    orthologs_dir = dir_vars["orthologs_parent"]
    formatted_df.to_csv(formatted_outpath, sep='\t')

    #Check for species from taxid_dict not represented in sql table
    other_tid_unique = orthologs_sql_df["Other_tax_id"].unique()
    absent_taxids = [taxid for taxid in taxid_dict if taxid not in other_tid_unique]
    return formatted_df, absent_taxids

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

def _compile_ortholog_data(local_orthologs_fpath, API_orthologs_fpath,tax_id_columns):
    local_orthologs = pd.read_csv(local_orthologs_fpath,sep='\t',dtype='str',index_col=0)
    API_orthologs = pd.read_csv(API_orthologs_fpath, sep='\t', dtype='str', index_col=0)
    local_missing_spec = ['29073']
    local_missing_cols = [tax_id_columns[col] for col in local_missing_spec]
    #Fill final_orthologs: 1) all local data, 2) add in API rows which aren't in local data
    # 3) correct empty cols (local_missing_spec) in local data with API data
    #3)
    union_index = local_orthologs.index.union(API_orthologs.index)
    final_orthologs = pd.DataFrame(index=union_index,columns=API_orthologs.columns)
    final_orthologs.loc[local_orthologs.index,:] = local_orthologs
    API_local_comp = API_orthologs.index.difference(local_orthologs.index)
    local_API_int = local_orthologs.index.intersection(API_orthologs.index)
    final_orthologs.loc[API_local_comp,:] = API_orthologs.loc[API_local_comp,:]
    final_orthologs.loc[local_API_int,local_missing_cols] = API_orthologs.loc[local_API_int,local_missing_cols]
    return final_orthologs

def _boolean_df_agg(ortholog_df,bool_df,how='any',comp=False):
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
    match_df = _boolean_df_agg(ortholog_table,bool_df,how=how,comp=comp)
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
    match_df = _boolean_df_agg(ortholog_table, bool_df, how=how, comp=comp)
    return match_df

def _identify_suspect_rows(local_orthologs_fpath,API_orthologs_fpath,match_pat,include_na=True):
    """

    :param local_orthologs_fpath:
    :param API_orthologs_fpath:
    :param include_na:
    :return: Rows from API orthologs not present in local orthologs which have pattern mismatches to match_pat
    """
    match_pat = "^10|^11"
    API_orthologs = load_orthologs_table(API_orthologs_fpath)
    print("API orthologs: {0}".format(len(API_orthologs)))
    local_orthologs = load_orthologs_table(local_orthologs_fpath)
    print("local orthologs: {0}".format(len(local_orthologs)))
    any_mismatch = ortholog_patmatch(API_orthologs,match_pat,how='all',comp=True)
    # API_all_valid = API_orthologs.dropna(how='any')
    if include_na:
        any_na = ortholog_na_filter(API_orthologs,how='any',comp=False)
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
    return _boolean_df_agg(sub_table,bool_df,how='any',comp=False)

def UCSC_NCBI_ID_query(driver,id_df,idx,hgid,tax_spec_names,tax_id_columns,
                       id_out_fpath="",initial_timeout=1):
    """For a human NCBI gene id, searches for ortholog gene ids corresponding to all taxids in taxid_columns.
    NOTE: Edited from spec_subs version of function to support multiple taxids per gene_id query

    :param driver (selenium.webdriver): Selenium webdriver object, initialized pre-iteration with headless option and implicit wait
    :param (DataFrame) id_df: DataFrame object where ortholog gene IDs will be stored
    :param (idx) idx: for row in id_df in which to store the ortholog gene ID
    :param (str) hgid: NCBI human Gene ID for given gene symbol
    :param (dict) tax_id_columns: dictionary from NCBI taxonomy IDs to column names in id_df
    :param (str) id_out_fpath: if provided, will write updated id_df as tsv to provided path.
    :param (int) initial_timeout: Initial timeout length for loading NCBI Gene orthologs page. Default 1

    :return: None. Edits id_df directly, or raises NCBIQueryError to log entries with no ortholog data for any species
    """

    #Mammalian taxonomy level: 40674 - ensure this is a valid taxonomy identifier and encompasses all NCBI species
    # in analysis or all requests will fail to produce ortholog data
    taxonomy_level = '40674'
    req_url = "https://www.ncbi.nlm.nih.gov/gene/{0}/ortholog/?scope={1}".format(hgid,taxonomy_level)
    driver.get(req_url)
    driver_timeout = initial_timeout

    for taxid in tax_spec_names:
        out_col = tax_id_columns[taxid]
        spec_name = tax_spec_names[taxid]
        gene_xpath = "//tr/td/label/*[contains(text(), '{0}')]".format(spec_name) + \
                     "/parent::label/parent::td/parent::tr/td[@class=' fld-gene']/a"
        try:
            spec_elems = WebDriverWait(driver,driver_timeout).until(
                EC.presence_of_all_elements_located((By.XPATH,gene_xpath))
            )
        except TimeoutException:
            #Check for generic ortholog table row elements (ie presence of any orthologs at all)
            generic_tr_xpath = "//tr/td/label"
            try:
                driver.find_element_by_xpath(generic_tr_xpath)
            except NoSuchElementException:
                # No orthologs at all in table, raise NCBIQueryError
                raise NCBIQueryError(0, "No orthologs available for NCBI hgid {0}".format(hgid))
            #If no exception raised, some orthologs exist, check for remaining
            id_df.loc[idx, out_col] = np.nan
            driver_timeout = 0.5
            continue
        if len(spec_elems) > 0:
            spec_href = spec_elems[0].get_attribute('href')
            spec_gid = re.search("/gene/(\d*)", spec_href).groups()[0]
            id_df.loc[idx,out_col] = spec_gid
            driver_timeout = 0
    if id_out_fpath:
        id_df.to_csv(id_out_fpath,sep='\t')

def map_AGS_geneIDs(xref_inpath, taxid_dict, results_outpath, errors_fpath,
                    overwrite_gid=[],tax_subset=None):
    """Using NCBI_gid entries storied in tsv_inpath, map ortholog IDs for species in taxid_dict
    from NCBI website and store in results_outpath.

    :param (str) xref_inpath: File path to xref table containing column "NCBI_gid" and indexded on UCSC_transcript_id
    :param (dictionary) taxid_dict: dict mapping NCBI Taxonomy IDs to species names. Species names must have
    conventional capitalization (ie 'Urocitellus parryii') and are used to match rows from NCBI query pages
    :param (str) results_outpath: File path for filled orthologs table.
    :param (str) errors_fpath: File path to errors log. Indexed on UCSC_tid for all project functions
    :param (array-like) overwrite_gid: If provided, will force queries for all hgid values in overwrite_gid
    :param (array-like) tax_subset: If provided, will only perform searches for taxids in tax_subset.
    Used in filling Ortholog table for species not present in local SQL gene_orthologs table
    :return: Updated version of out_tsv_df after all relevant NCBI Gene Ortholog queries
    """
    """Reads tsv file containing column of NCBI gene IDs, finds NCBI Gene IDs for all species in taxid_list
    """
    driver = headless_driver()
    xref_gid_col = "NCBI_gid"
    tax_gid_columns = ["{0}_gid".format(taxid) for taxid in taxid_dict]
    gid_columns_dict = dict(zip(taxid_dict.keys(),tax_gid_columns))
    if tax_subset:
        taxid_dict = {k:taxid_dict[k] for k in taxid_dict if k in tax_subset}
        tax_gid_columns = {k:tax_gid_columns[k] for k in tax_gid_columns if k in tax_subset}

    xref_df = load_NCBI_xref_table(xref_inpath,gid_dtype='str')

    #Load exisxting ortholog table or start new table if none exists at results_outpath
    if os.path.exists(results_outpath):
        out_tsv_df = pd.read_csv(results_outpath, dtype=str, delimiter='\t', index_col=0)
        check_existing = True
    else:
        out_tsv_df = pd.DataFrame(columns=tax_gid_columns)
        check_existing = False
    out_tsv_df.index.name = "NCBI_hgid"
    #Load existing errors_fpath
    check_errors, ncbi_errors_df = load_errors(errors_fpath,"NCBIQueryError")
    try:
        for ucsc_tid, row in xref_df.iterrows():
            NCBI_hgid = row[xref_gid_col]
            if out_tsv_df.index.dtype == 'str':
                idx = str(NCBI_hgid)
            elif out_tsv_df.index.dtype == 'int':
                idx = int(NCBI_hgid)
            """
             Prevent duplicate req_url attempts; Only perform query if 
             1) No out_tsv_df was read from file or NCBI_gid is specified in overwrite_gid
             2) out_tsv_df was read but NCBI_gid not in out_tsv_df
            """
            if (not check_existing or idx in overwrite_gid) or \
                (check_existing and idx not in out_tsv_df.index):
                #Skip query if previously caused NCBIQueryError (ie no orthologs listed)
                if check_errors and ucsc_tid in ncbi_errors_df["gene"].unique():
                    continue
                try:
                    UCSC_NCBI_ID_query(driver,out_tsv_df,idx,NCBI_hgid,taxid_dict,
                                   gid_columns_dict,id_out_fpath=results_outpath,initial_timeout=1)
                except NCBIQueryError as ncbi_e:
                    write_errors(errors_fpath,ucsc_tid,ncbi_e)

    except Exception as e:
        print(e)
        print("Stopping prematurely.")
        return out_tsv_df
    return out_tsv_df


def fetch_NCBI_protein_records(NCBI_gid_df, taxid_dict, dir_vars,
                               overwrite_fasta=[], NCBI_api_key='',tax_subset=None):
    """Downloads NCBI protein records for each NCBI Gene ID listed in AGS_gene_id_df.

    AGS_gene_id_df: DataFrame object with required columns 'Gene Symbol' and 'AGS Gene ID.' Gene symbol
    entries are used to name fasta files downloaded; AGS Gene IDs are queried using Entrez elink
    to 1) match Gene ID to all corresponding Protein IDs and 2) download those Protein IDs into one fasta
    file, saved into NCBI_homology_dirpath
    """

    orthologs_dir, orthologs_seq = dir_vars["orthologs_parent"],dir_vars["orthologs_seq"]
    try:
        ENTREZ_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        if NCBI_api_key:
            api_key_ext = "&api_key={0}".format(NCBI_api_key)
        else:
            api_key_ext = ''

        NCBI_protein_ID_tsv_fpath = "{0}/NCBI_protein_IDs.tsv".format(orthologs_dir)
        NCBI_protein_ID_columns = []
        for taxid in taxid_dict:
            gid_col = "{0}_gid".format(taxid)
            pid_col = "{0}_pids".format(taxid)
            NCBI_protein_ID_columns.extend([gid_col, pid_col])

        if tax_subset:
            taxid_dict = {k:taxid_dict[k] for k in taxid_dict if k in tax_subset}

        if os.path.exists(NCBI_protein_ID_tsv_fpath):
            NCBI_protein_ID_table = pd.read_csv(NCBI_protein_ID_tsv_fpath, dtype='str', delimiter='\t', index_col=0)
            NCBI_protein_ID_table.index.name = "NCBI_hgid"
        else:
            NCBI_protein_ID_table = pd.DataFrame(index=NCBI_gid_df.index, columns=NCBI_protein_ID_columns, dtype=str)
            NCBI_protein_ID_table.index.name = "NCBI_hgid"

            #     OVERWRITE_FASTAS = ["ATPIF1"]
        # Convert Gene ID to list of Protein IDs corresponding to transcript variant sequences
        for NCBI_gid, row in NCBI_gid_df.iterrows():
            # pid_row = NCBI_protein_ID_table.loc[NCBI_gid, :]
            for taxid in taxid_dict:
                gid_col = "{0}_gid".format(taxid)
                pid_col = "{0}_pids".format(taxid)
                spec_gid = row[gid_col]
                spec_dirpath = "{0}/{1}".format(orthologs_seq, taxid)
                spec_fasta_fpath = "{0}/{1}.fasta".format(spec_dirpath, NCBI_gid)
                #Check for initiating NCBI elink/efetch calls: If 1) spec_gid is not null and 2)
                #Fastas fetch not done previously or is being forced by overwrite_fasta
                if NCBI_gid in overwrite_fasta and os.path.exists(spec_fasta_fpath):
                    os.remove(spec_fasta_fpath)
                if (not os.path.exists(spec_fasta_fpath) or NCBI_gid in overwrite_fasta) \
                        and type(spec_gid) == str:
                    elink_req = "elink.fcgi?dbfrom=gene&db=protein&id={0}{1}".format(spec_gid, api_key_ext)
                    gp_elink_url = ENTREZ_BASE_URL + elink_req
                    elink_response = requests.get(gp_elink_url)
                    xml_data = elink_response.content
                    root = ET.fromstring(xml_data)
                    # Check XML formatting of elink pages - update xpath accordingly if functionality breaks
                    # Pulls Record IDs for Protein specifically; use gene_protein_refseq for Protein RefSeqs
                    protein_IDs = [link.text for link in root.findall(".//LinkSetDb[LinkName='gene_protein']/Link/Id")]
                    id_str = ','.join(protein_IDs)
                    NCBI_protein_ID_table.loc[NCBI_gid,gid_col] = spec_gid
                    NCBI_protein_ID_table.loc[NCBI_gid, pid_col] = id_str
                    efetch_req = "efetch.fcgi?db=protein&id={0}&rettype=fasta&retmode=text{1}".format(id_str,api_key_ext)
                    efetch_url = ENTREZ_BASE_URL + efetch_req
                    subprocess.run(args=['wget', efetch_url, '-O', spec_fasta_fpath])
                    NCBI_protein_ID_table.to_csv(NCBI_protein_ID_tsv_fpath, sep='\t')
        return NCBI_protein_ID_table
    except Exception as e:
        print(e)
        print("Updating to file and stopping prematurely")
        NCBI_protein_ID_table.to_csv(NCBI_protein_ID_tsv_fpath, sep='\t')
        return NCBI_protein_ID_table


def allseq_NCBI_UCSC_slignment(NCBI_xref_df, taxid_dict, dir_vars,gid_subset=[]):
    """File preparation for allNCBI data. Prepares: 1) compiled NCBI raw and aligned data and 2) all NCBI data profiled
    aligned against corresponding UCSC data.

    :param NCBI_xref_df:
    :param taxid_dict:
    :param dir_vars:
    :param gid_subset:
    :return:
    """
    from fastaUtility import MSA_fpath_list, profile_MSA
    # Directory paths
    dir_labels = ["orthologs_parent","orthologs_seq","UCSC_raw_parent","allNCBI_parent"]
    orthologs_dir, orthologs_seq, ucsc_raw, allNCBI_parent = [dir_vars[label] for label in dir_labels]
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
            NCBI_fpath = "{0}/{1}/{2}.fasta".format(orthologs_seq, taxid, NCBI_gid)
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
                    profile_MSA(in1_fpath=UCSC_raw_tid_fasta, in2_fpath=NCBI_alignment_fpath, \
                                out_fpath=combined_alignment_fpath)
                except AssertionError as e:
                    print("UCSC_tid_fasta_path: {0}".format(UCSC_raw_tid_fasta))
                    print("NCBI_alignment_fpath: {0}".format(NCBI_alignment_fpath))

        else:
            pass

def main():
    import directoryUtility
    from IPython.display import display
    config,taxid_dict,dir_vars = directoryUtility.config_initialization()
    orthologs_dir = dir_vars["orthologs_parent"]
    reorg_dir, dataset_name = dir_vars["reorg_parent"],dir_vars["dataset_name"]

    formatted_local_fpath = "{0}/knownCanonical_orthologs_local.tsv".format(orthologs_dir)
    local_orthologs_fpath = formatted_local_fpath
    API_orthologs_fpath = "{0}/NCBI_orthologs_API.tsv".format(orthologs_dir)
    NCBI_xref_inpath = "{0}/xref_summary_{1}/NCBI_xrefs.tsv".format(reorg_dir,dataset_name)
    dataset_id = dir_vars["dataset_identifier"]
    final_errors_fpath = "summary/{0}/errors.tsv".format(dataset_id)
    final_orthologs_fpath = "{0}/orthologs_final.tsv".format(orthologs_dir)
    sus_fpath = "{0}/NCBI_orthologs_031320.tsv".format(orthologs_dir)

    format_sql_check = False
    if format_sql_check:
        # format orthologs table
        formatted_df, absent_taxids = format_orthologs_sql(dir_vars, taxid_dict, formatted_local_fpath)
        display(formatted_df)

    id_mismatch_gids = False
    if id_mismatch_gids:
        valid_match_pat = "^10|^11"
        sus_ortho_df = load_orthologs_table()
        mismatch_df = ortholog_patmatch(sus_ortho_df, valid_match_pat, how='all', comp=True)
        with pd.option_context('display.max_columns', None):
            display(mismatch_df)
        suspect_entries_fpath = "{0}/NCBI_orthologs_suspect.tsv".format(orthologs_dir)
        mismatch_df.to_csv(suspect_entries_fpath, sep='\t')

    compile_final_orthologs = False
    if compile_final_orthologs:
        taxid_columns = dict(zip(taxid_dict.keys(),['{0}_gid'.format(k) for k in taxid_dict.keys()]))
        final_orthologs = _compile_ortholog_data(local_orthologs_fpath,API_orthologs_fpath,taxid_columns)
        final_orthologs.to_csv(final_orthologs_fpath,sep='\t')

    check_diff = True
    if check_diff:
        # diff_df = _compare_formatted_orthologs(sus_fpath,API_orthologs_fpath)
        diff_df = _compare_formatted_orthologs(sus_fpath,final_orthologs_fpath)
        display(diff_df)
        diff_df.to_csv('tmp/diff.tsv',sep='\t')

    # Optional repeat geneID mapping to fix na entries or just mismatch entries (ie invalid GeneID values).
    repeat_na_ID_mapping = False
    if repeat_na_ID_mapping:
        non_local_API_suspect = _identify_suspect_rows(local_orthologs_fpath,API_orthologs_fpath,valid_match_pat,include_na=True)
        map_AGS_geneIDs(NCBI_xref_inpath,taxid_dict,API_orthologs_fpath,final_errors_fpath,overwrite_gid=non_local_API_suspect.index)

    repeat_mismatch_mapping = False
    if repeat_mismatch_mapping:
        final_any_mismatch = ortholog_patmatch(final_orthologs,valid_match_pat,how='all',comp=True)
        map_AGS_geneIDs(NCBI_xref_inpath,taxid_dict,API_orthologs_fpath,final_errors_fpath,overwrite_gid=final_any_mismatch.index)

    # Update protein records for diff_df
    update_pids = False
    if update_pids:
        api_key = config["DIRECTORY"]["NCBIAPIKey"]
        pid_df = fetch_NCBI_protein_records(final_orthologs,taxid_dict,dir_vars,
                                        overwrite_fasta=diff_df.index,NCBI_api_key=api_key)

    process_allNCBI = True
    if process_allNCBI:
        NCBI_xrefs = load_NCBI_xref_table(NCBI_xref_inpath, gid_dtype='int')
        allseq_NCBI_UCSC_slignment(NCBI_xrefs,taxid_dict,dir_vars,gid_subset=diff_df.index.values)


if __name__ == "__main__":
    main()