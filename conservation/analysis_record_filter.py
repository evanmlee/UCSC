import pandas as pd
import numpy as np
import os
from IPython.display import display
from utility import fastaUtility as fastautil

def drop_NCBI_species_from_UCSC_records(combined_records_df,ncbi_index,analysis_taxid_dict):
    """

    :param combined_records_df: Records DataFrame containing UCSC record data and NCBI record(s)
    :param ncbi_index: Corresponds to NCBI species being used for analysis
    :param analysis_taxid_dict: Dictionary with NCBI taxids of analysis species as keys. UCSC records with any taxids
    in this set will be excluded from the UCSC used as the outgroup in analysis
    :return: filtered: combined_records_df with UCSC records corresponding to analysis_taxids removed
    """
    ncbi_partition = (combined_records_df.index.isin(ncbi_index)) & ~(combined_records_df.index.str.contains("ENST"))
    filtered_ucsc, ncbi_df = combined_records_df.loc[~ncbi_partition,:],combined_records_df.loc[ncbi_partition,:]
    for drop_taxid in analysis_taxid_dict:
        filtered_ucsc = filtered_ucsc.loc[filtered_ucsc['NCBI_taxid']!=drop_taxid,:]
    filtered = filtered_ucsc.append(ncbi_df)
    return filtered

def drop_redundant_UCSC_records(combined_records_df,ncbi_index):
    """For a given NCBI taxid, exclude UCSC from uniqueness/ conservation analysis from species which are either
    identical or close evolutionary relatives that likely contain similar adaptive residues to ncbi_taxid species
    (i.e. 13LGS and AGS). If provided NCBI taxids ar not in DROPPED_TAXIDS, returns unfilte
    :param combined_records_df: Records DataFrame containing UCSC record data and NCBI record(s)
    :param ncbi_index: Index or index value that corresponds to NCBI species for analysis.
    """
    DROPPED_TAXIDS = {9999:[43179],10181:[10181]}
    ncbi_taxids = combined_records_df.loc[ncbi_index,'NCBI_taxid']
    ncbi_partition = (combined_records_df.index.isin(ncbi_index)) & ~(combined_records_df.index.str.contains("ENST"))
    filtered_ucsc, ncbi_df = combined_records_df.loc[~ncbi_partition, :], combined_records_df.loc[ncbi_partition, :]
    for ncbi_taxid in ncbi_taxids:
        if ncbi_taxid in DROPPED_TAXIDS:
            dropped_set = DROPPED_TAXIDS[ncbi_taxid]
            filtered_ucsc = filtered_ucsc.loc[~filtered_ucsc['NCBI_taxid'].isin(dropped_set),:]
    filtered = filtered_ucsc.append(ncbi_df)
    return filtered

def filter_analysis_subset(combined_fasta,records_tsv_fpath,UCSC_analysis_subset=[],NCBI_record_subset=[],
                           filtered_outpath="tmp/filtered_analysis_set.fasta",taxid_dict={},
                           drop_redundant=True,drop_ncbi_from_ucsc=True):
    """Given subsets of UCSC and NCBI records to include in conservation_analysis, writes filtered sequence data to
    filtered_outpath and returns a DataFrame with filtered record information and the path to the sequence file.

    :param combined_fasta: bestNCBI alignment containing both all UCSC data and available selected NCBI records
    :param records_tsv_fpath: File path to for which analysis record subset information will be written
    :param UCSC_analysis_subset: If provided, limits UCSC records to those ids in UCSC_record_subset
    :param NCBI_record_subset: If provided, limits NCBI records to those ids in NCBI_record_subset
    :param filtered_outpath: optional parameter. If provided, writes filtered fasta to this path. Otherwise writes to
    a tmp record file.
    :param taxid_dict: If provided, will be used instead of config/NCBI_analysis_tax.txt when reading NCBI
    record information.
    :param drop_redundant: Default True. Removes UCSC records from final analysis set which correspond to identical/
    extremely close evo relatives for sake of determining unique residues/ analysis for NCBI species.
    :return: records_df: DataFrame containing record information represented in filtered sequence set
    :return: filtered_outpath: Outpath to which filtered records were written
    """

    records_df = fastautil.load_UCSC_NCBI_df(combined_fasta,ncbi_taxid_dict=taxid_dict,
                                         UCSC_subset=UCSC_analysis_subset,NCBI_subset=NCBI_record_subset)
    if drop_redundant:
        records_df = drop_redundant_UCSC_records(records_df,ncbi_index=NCBI_record_subset)
    if drop_ncbi_from_ucsc:
        records_df = drop_NCBI_species_from_UCSC_records(records_df,NCBI_record_subset,taxid_dict)
    records_df.drop(columns=['sequence'],inplace=True)
    records_df.to_csv(records_tsv_fpath,sep='\t')
    fastautil.filter_fasta_infile(records_df.index,combined_fasta,outfile_path=filtered_outpath)
    return records_df, filtered_outpath

def length_metric_supplement(taxid_dict,dir_vars,length_metrics_fpath,force_overwrite=False,
                             suppl_lm_df_fpath=""):
    """Adds data to length_metrics table corresponding to record lengths for species in UCSC_taxids and NCBI_taxids.
    Writes updated table to suppl_lm_df_fpath.

    :param taxid_dict: Maps NCBI taxids to species names for species in NCBI analysis set
    :param dir_vars: Contains key 'bestNCBI_parent' and value corresponding to path to bestNCBI directory
    :param length_metrics_fpath: pre-existing length-metrics table file path
    :param force_overwrite: If true, will overwrite file at suppl_lm_df_fpath even if it exists
    :param suppl_lm_df_fpath: If not provided, defaults to length_metrics_suppl.tsv in same directory as
    length_metrics_fpath
    :return: suppl_lm_df_fpath: File path to supplemented length metrics table
    """
    from IPython.display import display
    lm_df = pd.read_csv(length_metrics_fpath,sep='\t',index_col=0)
    if not suppl_lm_df_fpath:
        suppl_lm_df_fpath = length_metrics_fpath[:-4]+"_suppl.tsv"
    bestNCBI_dir = dir_vars['bestNCBI_parent']
    combined_dir = "{0}/combined".format(bestNCBI_dir)
    UCSC_taxids = [9606,10090,43179,10181]
    UCSC_cols = ['9606_length','10090_length','43179_length','10181_UCSC_length']
    NCBI_taxids = [9999,10181,29073,9994]
    NCBI_cols = ['9999_length','10181_length','29073_length','9994_length']
    if os.path.exists(suppl_lm_df_fpath) and not force_overwrite:
        print("File already exists at supplemented length metrics path: {0}".format(suppl_lm_df_fpath))
        print("Run with force_overwrite=True if recalculation desired.")
        return suppl_lm_df_fpath
    for tid,row in lm_df.iterrows():
        comb_fasta = "{0}/{1}.fasta".format(combined_dir,tid)
        comb_df = fastautil.load_UCSC_NCBI_df(comb_fasta,taxid_dict)
        for i,taxid in enumerate(UCSC_taxids):
            spec_length = comb_df.loc[(comb_df['NCBI_taxid']==taxid) & (comb_df.index.str.contains("ENST")), 'length']
            if len(spec_length) == 1:
                lm_df.loc[tid,UCSC_cols[i]] = spec_length.iloc[0]
            elif len(spec_length) == 0:
                continue
            else:
                print("More than one length value for taxid.")
                display(spec_length)
                print("Transcript ID: {0}".format(tid))
                print("TaxID: {0}".format(taxid))
        for i,taxid in enumerate(NCBI_taxids):
            spec_length = comb_df.loc[(comb_df['NCBI_taxid']==taxid) & ~(comb_df.index.str.contains("ENST")), 'length']
            if len(spec_length) == 1:
                lm_df.loc[tid,NCBI_cols[i]] = spec_length.iloc[0]
            elif len(spec_length) == 0:
                continue
            else:
                print("More than one length value for taxid.")
                display(spec_length)
                print("Transcript ID: {0}".format(tid))
                print("TaxID: {0}".format(taxid))
    lm_df.to_csv(suppl_lm_df_fpath,sep='\t')
    return suppl_lm_df_fpath

def length_metric_filter(taxid_dict,dir_vars,length_metrics_fpath,tolerance=0.05,how='all'):
    """Returns a DataFrame check_df with transcript_id index and columns corresponding to taxonomy IDs and values
    on whether that species record passed length checks.

    :param (dict) taxid_dict: maps taxid to species name for NCBI analysis species
    :param dir_vars: maps 'bestNCBI_parent' to directory path for bestNCBI-record-against-UCSC alignment data
    :param length_metrics_fpath: file path to length_metrics_table
    :param tolerance: Maximum percent difference in length tolerated for a sequence record to pass length checks
    :param (str) how: 'any' or 'all'
    :return:
    Length checks are against 1) homo sapiens record, 2)
    """
    suppl_lm_df_fpath = length_metric_supplement(taxid_dict,dir_vars,length_metrics_fpath)
    check_col_labels = ["{0}_check".format(taxid) for taxid in taxid_dict]
    check_df = pd.DataFrame(columns=check_col_labels)
    suppl_lm_df = pd.read_csv(suppl_lm_df_fpath,sep='\t',index_col=0)
    suppl_row_generator = suppl_lm_df.iterrows()
    tax_specific_checks = {9999:['43179_length'],10181:['10181_UCSC_length'],29073:[],9994:[]}

    for taxid in taxid_dict:
        length_label = "{0}_length".format(taxid)
        check_col = "{0}_check".format(taxid)
        check_labels = ['euth_median','bestncbi_median','9606_length','10090_length']
        check_labels.extend(tax_specific_checks[taxid])
        idx, row = suppl_row_generator.__next__()
        spec_indiv_checks = pd.DataFrame(columns=check_labels)
        # print(idx)
        # for idx, row in suppl_row_generator:
        if True:
            spec_length = row[length_label]
            if np.isnan(spec_length):
                check_df.loc[idx,check_col] = np.nan#False#np.nan #False
                spec_indiv_checks.loc[idx,:] = np.nan
                continue
            check_lengths = row[check_labels]
            checks = (check_lengths-spec_length)/check_lengths <= tolerance
            if how == 'any':
                check_df.loc[idx,check_col] = checks.any()
                spec_indiv_checks.loc[idx, :] = checks
            elif how=='all':
                check_df.loc[idx, check_col] = checks.all()
                spec_indiv_checks.loc[idx, :] = checks

    display(check_df)
