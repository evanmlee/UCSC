import pandas as pd
import os

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