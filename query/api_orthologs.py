import pandas as pd
import os
import re
import numpy as np
import subprocess
import requests
import xml.etree.ElementTree as ET

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException,NoSuchElementException
from utility.UCSCerrors import NCBIQueryError, write_errors, load_errors

from query.orthologUtility import load_NCBI_xref_table

ENTREZ_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def headless_driver():
    #Headless driver with default options
    from selenium.webdriver.chrome.options import Options
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    WINDOW_SIZE = "1920,1080"
    chrome_options.add_argument("--window-size=%s" % WINDOW_SIZE)
    driver = webdriver.Chrome(options=chrome_options)
    return driver

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
                if check_errors and ucsc_tid in ncbi_errors_df["tid"].unique():
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

def species_ortholog_request(NCBI_protein_ID_table,spec_fasta_fpath,spec_gid,taxid,NCBI_hgid,api_key_ext='',
                             pid_tsv_fpath=""):
    gid_col = "{0}_gid".format(taxid)
    pid_col = "{0}_pids".format(taxid)
    elink_req = "elink.fcgi?dbfrom=gene&db=protein&id={0}{1}".format(spec_gid, api_key_ext)
    gp_elink_url = ENTREZ_BASE_URL + elink_req
    elink_response = requests.get(gp_elink_url)
    xml_data = elink_response.content
    root = ET.fromstring(xml_data)
    # Check XML formatting of elink pages - update xpath accordingly if functionality breaks
    # Pulls Record IDs for Protein specifically; use gene_protein_refseq for Protein RefSeqs
    protein_IDs = [link.text for link in root.findall(".//LinkSetDb[LinkName='gene_protein']/Link/Id")]
    id_str = ','.join(protein_IDs)
    NCBI_protein_ID_table.loc[NCBI_hgid, gid_col] = spec_gid
    NCBI_protein_ID_table.loc[NCBI_hgid, pid_col] = id_str
    efetch_req = "efetch.fcgi?db=protein&id={0}&rettype=fasta&retmode=text{1}".format(id_str, api_key_ext)
    efetch_url = ENTREZ_BASE_URL + efetch_req
    subprocess.run(args=['wget', efetch_url, '-O', spec_fasta_fpath])
    if pid_tsv_fpath:
        NCBI_protein_ID_table.to_csv(pid_tsv_fpath, sep='\t')
    return NCBI_protein_ID_table

def fetch_NCBI_protein_records(NCBI_gid_df, taxid_dict, dir_vars,
                               overwrite_fasta=[], NCBI_api_key='',tax_subset=None):
    """Downloads NCBI protein records for each NCBI Gene ID listed in AGS_gene_id_df.

    AGS_gene_id_df: DataFrame object with required columns 'Gene Symbol' and 'AGS Gene ID.' Gene symbol
    entries are used to name fasta files downloaded; AGS Gene IDs are queried using Entrez elink
    to 1) match Gene ID to all corresponding Protein IDs and 2) download those Protein IDs into one fasta
    file, saved into NCBI_homology_dirpath
    """

    orthologs_dir, orthologs_seq = dir_vars["orthologs_parent"],dir_vars["orthologs_seq"]
    if NCBI_api_key:
        api_key_ext = "&api_key={0}".format(NCBI_api_key)
    else:
        api_key_ext = ''

    protein_ID_tsv_fpath = "{0}/NCBI_protein_IDs.tsv".format(orthologs_dir)
    NCBI_protein_ID_columns = []
    for taxid in taxid_dict:
        gid_col = "{0}_gid".format(taxid)
        pid_col = "{0}_pids".format(taxid)
        NCBI_protein_ID_columns.extend([gid_col, pid_col])
    if tax_subset:
        taxid_dict = {k:taxid_dict[k] for k in taxid_dict if k in tax_subset}
    if os.path.exists(protein_ID_tsv_fpath):
        NCBI_protein_ID_table = pd.read_csv(protein_ID_tsv_fpath, dtype='str', delimiter='\t', index_col=0)
    else:
        NCBI_protein_ID_table = pd.DataFrame(index=NCBI_gid_df.index, columns=NCBI_protein_ID_columns, dtype=str)
    NCBI_protein_ID_table.index.name = "NCBI_hgid"

    # Convert Gene ID to list of Protein IDs corresponding to transcript variant sequences
    for NCBI_hgid, row in NCBI_gid_df.iterrows():
        # pid_row = NCBI_protein_ID_table.loc[NCBI_gid, :]
        for taxid in taxid_dict:
            gid_col = "{0}_gid".format(taxid)
            spec_gid = row[gid_col]
            spec_fasta_fpath = "{0}/{1}/{2}.fasta".format(orthologs_seq,taxid,NCBI_hgid)
            #Check for initiating NCBI elink/efetch calls: If 1) spec_gid is not null and 2)
            #Fastas fetch not done previously or is being forced by overwrite_fasta
            if NCBI_hgid in overwrite_fasta and os.path.exists(spec_fasta_fpath):
                os.remove(spec_fasta_fpath)
            if (not os.path.exists(spec_fasta_fpath) or NCBI_hgid in overwrite_fasta) \
                    and type(spec_gid) == str:
                species_ortholog_request(NCBI_protein_ID_table,spec_fasta_fpath,spec_gid,taxid,NCBI_hgid,
                                         api_key_ext=api_key_ext,pid_tsv_fpath=protein_ID_tsv_fpath)
