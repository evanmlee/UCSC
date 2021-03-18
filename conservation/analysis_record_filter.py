import pandas as pd
import numpy as np
import os
import warnings
from IPython.display import display

from utility import fastaUtility as fautil, UCSCerrors
from utility import NCBIfilter as filt
from utility.directoryUtility import taxid_dict,dir_vars,ucsc_taxid_dict
from utility.NCBIfilter import CLOSEST_EVO_TAXIDS, SECONDARY_EVO_TAXIDS

def drop_NCBI_species_from_UCSC_records(filtered_ucsc,analysis_taxid_dict):
    """Drops UCSC record data for any species tax_ids represented in analysis_taxid_dict.
    :param combined_records_df: Records DataFrame containing UCSC record data and NCBI record(s)
    :param ncbi_idx: Corresponds to NCBI species being used for analysis
    :param analysis_taxid_dict: Dictionary with NCBI taxids of analysis species as keys. UCSC records with any taxids
    in this set will be excluded from the UCSC used as the outgroup in analysis
    :return: filtered: combined_records_df with UCSC records corresponding to analysis_taxids removed
    """
    for drop_taxid in analysis_taxid_dict:
        filtered_ucsc = filtered_ucsc.loc[filtered_ucsc['NCBI_taxid']!=drop_taxid,:]
    return filtered_ucsc

def drop_ucsc_analysis_species(filtered_ucsc):
    """ Drops UCSC data for species in [UCSC_TAXONOMY] section of config.txt from the mammalian outgroup DataFrame.
    :param filtered_ucsc: Mammalian outgroup DataFrame (test record should have been removed in filter_analysis_subset)
    :return: filtered_ucsc: corresponding DataFrame with UCSC analysis species removed
    """
    for drop_taxid in ucsc_taxid_dict:
        filtered_ucsc = filtered_ucsc.loc[filtered_ucsc['NCBI_taxid'] != drop_taxid, :]
    return filtered_ucsc

def drop_incomplete_records(ucsc_df, test_df,length_thresh=0.7, how='test_record'):
    """Applies a length filter to combined_records_df and drops UCSC records which have length less than accepted record
    lengths multiplied by length_thresh.

    :param combined_records_df: DataFrame containing record data from UCSC/NCBI. 
    :param ncbi_idx: Index of NCBI records in combined_records_df. Assumed to only have one record.
    :param (float) length_thresh: Should be between 0 and 1.
    :param how: 'ncbi': Uses records in NCBI idx as accepted records to set length threshold
                'nonzero_ucsc': Uses non-zero length UCSC records to determine median for length threshold
                'nonzero_evo': Uses non-zero evolutionary relative. If none, auto defaults to ncbi record length
    :return:
    """
    accepted_how_values = ['test_record','nonzero_ucsc','nonzero_evo']
    if how == 'test_record':
        accepted_records = test_df.index
        accepted_med_len = np.median(test_df.loc[accepted_records, 'length'])
    elif how == 'nonzero_ucsc':
        non_zero_ucsc = ucsc_df.loc[ucsc_df['length'] > 0, :]
        accepted_records = non_zero_ucsc.index
        accepted_med_len = np.median(ucsc_df.loc[accepted_records, 'length'])
    elif how == 'nonzero_evo':
        ncbi_taxid = test_df["NCBI_taxid"].iloc[0]
        evo_taxids = [CLOSEST_EVO_TAXIDS[ncbi_taxid]] + SECONDARY_EVO_TAXIDS[ncbi_taxid]
        evo_ucsc = ucsc_df.loc[ucsc_df['NCBI_taxid'].isin(evo_taxids),:]
        nonzero_evo = evo_ucsc.loc[evo_ucsc['length']> 0,:]
        if len(nonzero_evo) > 0:
            accepted_records = nonzero_evo.index
            accepted_med_len = np.median(ucsc_df.loc[accepted_records, 'length'])
        else:
            accepted_records = test_df.index
            accepted_med_len = np.median(test_df.loc[accepted_records, 'length'])
    else:
        raise ValueError("how must use values in accepted values: {0}".format((accepted_how_values)))
    cutoff_length = accepted_med_len*length_thresh
    filtered_ucsc_df = ucsc_df.loc[ucsc_df['length']>=cutoff_length,:]
    return filtered_ucsc_df

def filter_ucsc_df(fasta_fpath,ucsc_df,test_df,drop_ucsc_analysis_specs=True,drop_ncbi_from_ucsc=True,drop_incomplete=True):
    """
    :param fasta_fpath:
    :param ucsc_df: Unfiltered mammalian outgroup df.
    :param test_df:
    :param drop_ucsc_analysis_specs: Default True. Drops UCSC data for species in [UCSC_TAXONOMY] section of config.txt
    from the mammalian outgroup
    :param drop_ncbi_from_ucsc: Default True. Drops UCSC data corresponding to any species in the NCBI analysis set
    (i.e. if heterocephalus glaber is in NCBI species set, corresponding UCSC hetGla records will be dropped from
    outgroup for all NCBI analysis species).
    :param drop_incomplete: Default True. Applies length filter to UCSC records and drops any record below length
    threshold (determined by 80% of NCBI record length by default)
    :return:
    """
    if drop_ucsc_analysis_specs:
        ucsc_df = drop_ucsc_analysis_species(ucsc_df)
    if drop_ncbi_from_ucsc:
        ucsc_df = drop_NCBI_species_from_UCSC_records(ucsc_df,taxid_dict)
    if drop_incomplete:
        ucsc_df = drop_incomplete_records(ucsc_df,test_df, how='nonzero_ucsc')
    if len(ucsc_df) == 0:
        emsg = "{0}".format("Analysis records filtering resulted in empty outgroup set.")
        raise UCSCerrors.SequenceAnalysisError(0,emsg)
    elif len(ucsc_df)< 5:
        wmsg = "Fewer than 5 records in outgroup dataset after filtering from {0}".format(fasta_fpath)
        warnings.warn(wmsg,RuntimeWarning)
    return ucsc_df
    
def filter_analysis_subset(fasta_fpath,records_tsv_fpath,test_taxid,test_source="UCSC",UCSC_analysis_subset=[],
                           NCBI_record_subset=[],
                           filtered_outpath="tmp/filtered_analysis_set.fasta",drop_ucsc_analysis_specs=True,
                           drop_ncbi_from_ucsc=True,drop_incomplete=True):
    """Given subsets of UCSC and NCBI records to include in conservation_analysis, writes filtered sequence data to
    filtered_outpath and returns a DataFrame with filtered record information and the path to the sequence file.

    :param fasta_fpath: bestNCBI alignment containing both all UCSC data and available selected NCBI records, or
            clade-filtered UCSC data if test_taxid corresponds to a UCSC data species
    :param records_tsv_fpath: File path to for which analysis record subset information will be written
    :param UCSC_analysis_subset: If provided, limits UCSC records to those ids in UCSC_record_subset
    :param NCBI_record_subset: If provided, limits NCBI records to those ids in NCBI_record_subset
    :param filtered_outpath: optional parameter. If provided, writes filtered fasta to this path. Otherwise writes to
    a tmp record file.
    :param taxid_dict: If provided, will be used instead of config/NCBI_analysis_tax.txt when reading NCBI
    record information.
    :param drop_redundant: Default True. Removes UCSC records from final analysis set which correspond to identical/
    extremely close evo relatives for sake of determining unique residues/ analysis for NCBI species.
    :param drop_ncbi_from_ucsc: Default True. Drops UCSC data corresponding to any species in the NCBI analysis set
    (i.e. if heterocephalus glaber is in NCBI species set, corresponding UCSC hetGla records will be dropped from
    outgroup for all NCBI analysis species).
    :param drop_incomplete: Default True. Applies length filter to UCSC records and drops any record below length
    threshold (determined by 80% of NCBI record length by default)
    :return: records_df: DataFrame containing record information represented in filtered sequence set
    :return: filtered_outpath: Outpath to which filtered records were written
    """
    if test_source == "UCSC":
        records_df = fautil.load_UCSC_fasta_df(fasta_fpath)
        test_partition = (records_df["NCBI_taxid"]==test_taxid)
    elif test_source == "NCBI":
        records_df = fautil.load_UCSC_NCBI_df(fasta_fpath,ncbi_taxid_dict=taxid_dict,
                                         UCSC_subset=UCSC_analysis_subset,NCBI_subset=NCBI_record_subset)
        test_partition = (~records_df.index.str.contains("ENST")) & (records_df["NCBI_taxid"]==test_taxid)
    else:
        raise ValueError("Specify test_source as 'UCSC' or 'NCBI'")
    outgroup_ucsc_df,test_df = records_df.loc[~test_partition],records_df.loc[test_partition,:]
    filtered_ucsc = filter_ucsc_df(fasta_fpath,outgroup_ucsc_df,test_df,
                                   drop_ucsc_analysis_specs,drop_ncbi_from_ucsc,drop_incomplete)
    records_df = filtered_ucsc.append(test_df)
    records_df.drop(columns=['sequence'],inplace=True)
    records_df.to_csv(records_tsv_fpath,sep='\t')
    if test_source == 'UCSC':
        fautil.filter_fasta_infile(records_df.index, fasta_fpath, outfile_path=filtered_outpath,ordered=True)
    else:
        fautil.filter_fasta_infile(records_df.index,fasta_fpath,outfile_path=filtered_outpath)
    return records_df, filtered_outpath

def id_length_check(combined_fasta_fpath,lm_row,check_taxids,len_tol=0.05,min_id_tol=0.15,max_id_tol=0.35,
                    ce_id_adjust=1.25,quiet=True):
    """For a given species represented by taxid, runs length - identity checks for record validity and returns a dict
    mapping check labels to True/ False values.

    Length alone iscompared to the clade median and bestncbi median present in lm_row. For each species taxid in
    check_taxids, checks that either length of ncbi record is within len_tol (5%) of check_taxid species length. If
    this is not true, compares record identity against the check_taxid species. A correction is made for the length
    discrepancy between the ucsc and ncbi records and if the corrected identity is below id_tol, then the check will
    pass.
    Any record for check_taxids species with a length of 0 (i.e. no sequence data) will not compute a boolean
    check_value and should not count toward the total number of record checks (see use in length_metrics_checks).
    :param combined_fasta_fpath: fasta path to bestNCBI/combined alignment (single best NCBI ortholog record against
    raw UCSC clade filtered data).
    :param lm_row: length_metrics.tsv equivalent DataFrame indexed on transcript IDs. See utility/NCBIfilter.py
    :param check_taxids: Collection of NCBI taxonomy IDs corresponding to the correct NCBI species from combined_fasta
    for which length/ identity checks will be computed.
    :param len_tol: Tolerance for difference in record length. For any record which length corrected difference <=
    tolerance, corresponding check will be True
    :param id_tol: Tolerance for length corrected identity. If a record has a length difference over len_tol but
    a length corrected identity value below id_tol, corresponding check will be True
    :return: spec_checks: Dict, maps check labels (ie 'clade_check', '43179_check') to True/ False values. Records
    checks for which no ucsc_record is present (length is 0) will not be in spec_checks and will be ignored for the sake
    of checking overall check pass ratio.
    """
    combined_df = fautil.load_UCSC_NCBI_df(combined_fasta_fpath,taxid_dict)
    ucsc_mask = combined_df.index.str.contains("ENST")
    ucsc_df,ncbi_df = combined_df.loc[ucsc_mask],combined_df.loc[~ucsc_mask]

    ncbi_len = ncbi_df.loc[:,"length"][0]
    #Check clade median and best_NCBI_median
    clade_med, bestncbi_med = lm_row['clade_median'],lm_row['bestncbi_median']
    clade_check= (clade_med > 0 and np.abs((clade_med-ncbi_len)/clade_med)<=len_tol)
    bestncbi_check = np.abs((bestncbi_med-ncbi_len)/bestncbi_med)<=len_tol
    spec_checks = {"clade_check":clade_check,"bestncbi_check":bestncbi_check}
    precalc_id_dm = False

    for check_taxid in check_taxids:
        length_label, check_label = "{0}_length".format(check_taxid),"{0}_check".format(check_taxid)
        check_len = lm_row[length_label]
        if check_len == 0:
            continue
        length_discrepancy = (check_len-ncbi_len)/check_len
        #If needed, calculate id_dm and associated variables once. Filter file to speed up id calculation time
        if not precalc_id_dm:
            filtered_ids = ucsc_df.loc[ucsc_df['NCBI_taxid'].isin(check_taxids)].index.append(ncbi_df.index)
            tmp_fpath = "tmp/lm_checks.fasta"
            fautil.filter_fasta_infile(filtered_ids, combined_fasta_fpath, tmp_fpath)
            id_dm, align_srs = fautil.construct_id_dm(filtered_ids, tmp_fpath, filtered=True, aligned=True)
            ncbi_pos = align_srs.index.get_loc(ncbi_df.index[0])
            ncbi_id_row = id_dm[ncbi_pos, :]
            closest_taxid = check_taxids[0]
            ce_record = ucsc_df.loc[ucsc_df['NCBI_taxid']==closest_taxid,:].index[0]
            ce_row = id_dm[align_srs.index.get_loc(ce_record),:]
            precalc_id_dm = True
            if not quiet:
                print("ID DM: {0}".format(id_dm))
                print("DM Record Labels: {0}".format(align_srs.index.values))

        check_record = ucsc_df.loc[ucsc_df['NCBI_taxid']==check_taxid].index[0]
        check_dm_pos = align_srs.index.get_loc(check_record)
        overall_id_val = ncbi_id_row[check_dm_pos]
        #Correct for terminal Z character in UCSC sequences
        align_len = len(align_srs.iloc[0])
        ucsc_id_len_correction = 1/align_len
        #Aligned_portion_identity: Subtract out length discrepancy fraction of longest record length from id_dm value
        #Scale by inverse fraction of non-discrepancy length (ie max alignable portion given sequence lengths)
        align_discr_ratio = np.abs(check_len-ncbi_len)/align_len
        aligned_check_ratio = (1 - align_discr_ratio)
        aligned_portion_identity = (overall_id_val - align_discr_ratio  - ucsc_id_len_correction)/aligned_check_ratio
        ce_check_id = ce_id_adjust*(ce_row[check_dm_pos]/aligned_check_ratio)
        id_tol = min(max_id_tol,max(ce_check_id,min_id_tol))
        if not quiet:
            print("Check TaxID: {0}".format(check_taxid))
            print("ID tolerance threshold: {0}".format(id_tol))
            print("Check TaxID Against Closest Evo ID: {0}".format(ce_row[check_dm_pos]))
            print("Length adjusted identity: {0}".format(aligned_portion_identity))
        spec_checks[check_label] = (aligned_portion_identity <= id_tol)
    return spec_checks

def write_lm_params_file(params):
    params_fpath = "{0}/length_checks_params.txt".format(dir_vars['combined_run'])
    with open(params_fpath,'wt') as params_f:
        for label in params:
            val = params[label]
            fline = "{0}: {1}\n".format(label,val)
            params_f.write(fline)

def length_checks_df(lc_fpath,taxids,check_cols=True):
    """Reads length checks DataFrame from file if exists, otherwise generates empty DataFrame with columns based
    on taxids."""
    check_col_labels = ["{0}_check".format(taxid) for taxid in taxids]
    if os.path.exists(lc_fpath):
        check_df = pd.read_csv(lc_fpath,sep='\t',index_col=0)
        if check_cols:
            for col in check_df.columns:
                if col not in check_col_labels:
                    raise RuntimeError("Column labels for taxids ({0}) do not match columns list at specified file {1}"\
                                       .format(taxids, lc_fpath))
            for col in check_col_labels:
                if col not in check_df.columns:
                    raise RuntimeError("Column labels for taxids ({0}) do not match columns list at specified file {1}"\
                                       .format(taxids, lc_fpath))
    else:
        check_df = pd.DataFrame(columns=check_col_labels)
    check_df.index.name = "UCSC_transcript_id"
    return check_df


def length_metric_checks(length_metrics_fpath,len_tol=0.05,min_id_tol=0.15,max_id_tol=0.35,ce_id_adjust=1.25,
                         evo_how='all',pass_how='most',
                         recalc_spec_checks=False,recalc_pass_checks=False):
    """Returns a DataFrame check_df with transcript_id index and columns corresponding to taxonomy IDs and values
    on whether that species record passed length checks.

    :param length_metrics_fpath: file path to length_metrics_table
    :param len_tol: Maximum percent difference in length tolerated for a sequence record to pass length checks
    :param min_id_tol, max_id_tol: If NCBI record length doesn't meet len_tol, check length-adjusted identity. Threshold
     identity value used to pass will be based on closest_evo record distance to check record. If no closest evo record
     is present
    :param (str) evo_how: accepted values: 'all' or 'closest'. Determines which evolutionary relatives are used to
    calculate species specific checks. 'all': closest and secondary, 'closest': closest only
    :param (str) pass_how: accepted values: 'any','all','most'. Any: any length checks pass->True. All: all checks->True.
    Most: Over half pass->True
    :param recalc_spec_checks:
    :return: lm_checks: DataFrame, for each transcript_id in length metrics table, containsboolean values for each
    taxonomy ID which correspond to whether the bestNCBI record for that species passed length checks as specified with
    how.
    Length checks are against 1) boreoeutherian UCSC median 2) best NCBI records median 3) homo sapiens record,
    4) mus musculus record and 5) For heterocephalus glaber: checks against UCSC data for het g. For Urocitellus parryii:
    checks against speTri record length in UCSC data.
    """
    closest_evo_taxids, secondary_evo_taxids = filt.CLOSEST_EVO_TAXIDS, filt.SECONDARY_EVO_TAXIDS
    evo_taxids = {}
    for taxid in taxid_dict:
        if evo_how == 'all':
            check_taxids = [closest_evo_taxids[taxid]] + secondary_evo_taxids[taxid]
        elif evo_how == 'closest':
            check_taxids = [closest_evo_taxids[taxid]]
        else:
            raise ValueError("parameter 'evo_how' must be 'all' or 'closest'.")
        evo_taxids[taxid] = check_taxids
    checks_fpath = "{0}/length_checks.tsv".format(dir_vars['combined_run'])
    if os.path.exists(checks_fpath) and not recalc_pass_checks:
        return
    else:
        check_df = length_checks_df(checks_fpath,taxid_dict)
    lm_df,_ = filt.length_metrics_df(length_metrics_fpath,taxid_dict)
    for taxid in taxid_dict:
        tax_length_label = "{0}_ncbi_length".format(taxid)
        check_col = "{0}_check".format(taxid)
        check_labels = ["clade_check", "bestncbi_check"]
        check_taxids = evo_taxids[taxid]
        tax_specific_check_labels = ["{0}_check".format(check_taxid) for check_taxid in check_taxids]
        check_labels.extend(tax_specific_check_labels)
        spec_indiv_checks_fpath = "{0}/{1}_length_checks.tsv".format(dir_vars['combined_run'], taxid)
        if os.path.exists(spec_indiv_checks_fpath):
            spec_indiv_checks = pd.read_csv(spec_indiv_checks_fpath,sep='\t',index_col=0)
        else:
            spec_indiv_checks = pd.DataFrame(columns=check_labels)
        for ucsc_tid,lm_row in lm_df.iterrows():
            spec_length = lm_row[tax_length_label]
            if np.isnan(spec_length):
                #If no spec_length (ie no available record for species) use np.nan as fill values to differentiate
                #from failed checks
                check_df.loc[ucsc_tid,check_col] = np.nan
                spec_indiv_checks.loc[ucsc_tid,:] = np.nan
            else:
                if ucsc_tid in spec_indiv_checks.index and not recalc_spec_checks:
                    spec_row = spec_indiv_checks.loc[ucsc_tid,:]
                    nonna = spec_row.dropna()
                    pass_ratio = sum(nonna)/len(nonna)
                    ncbi_only_check = (nonna['bestncbi_check'] and sum(nonna) == 1)
                else:
                    combined_fasta_fpath = "{0}/combined/{1}/{2}.fasta".format(dir_vars['bestNCBI_parent'],taxid,ucsc_tid)
                    spec_checks = id_length_check(combined_fasta_fpath,lm_row,check_taxids,len_tol,min_id_tol,max_id_tol,
                                                  ce_id_adjust=ce_id_adjust)
                    spec_indiv_checks.loc[ucsc_tid, :] = spec_checks
                    pass_ratio = sum(spec_checks.values())/len(spec_checks)
                    ncbi_only_check = (spec_checks['bestncbi_check']==True and sum(spec_checks.values()) == 1)
                if pass_how == 'any':
                    pass_check = (pass_ratio > 0) and not ncbi_only_check
                    check_df.loc[ucsc_tid,check_col] = pass_check
                elif pass_how=='all':
                    pass_check = (pass_ratio == 1) and not ncbi_only_check
                    check_df.loc[ucsc_tid, check_col] = pass_check
                elif pass_how=='most':
                    pass_check = (pass_ratio >= 0.5) and not ncbi_only_check
                    check_df.loc[ucsc_tid, check_col] = pass_check
        if recalc_spec_checks:
            spec_indiv_checks.to_csv(spec_indiv_checks_fpath,sep='\t')

    check_df.to_csv(checks_fpath,sep='\t')
    param_labels = ["min_id_tol","max_id_tol","closest_evo_thresh_adjust","evo_how","pass_how"]
    param_vals = [min_id_tol,max_id_tol,ce_id_adjust,evo_how,pass_how]
    params = dict(zip(param_labels,param_vals))
    write_lm_params_file(params)


def ucsc_length_checks(clade_level="boreoeutheria",len_tol=0.7,length_how="non_zero_clade",
                       recalc_pass_checks=False,alt_ucsc_lc_fpath=""):
    import query.orthologUtility as orutil
    any_ncbi_ids = orutil.any_NCBI_xref()
    reorg_parent, ds_id = dir_vars['reorg_parent'],dir_vars['dataset_identifier']
    if alt_ucsc_lc_fpath:
        ucsc_lc_fpath = alt_ucsc_lc_fpath
    else:
        ucsc_lc_fpath = "{0}/ucsc_length_checks.tsv".format(dir_vars['combined_run'])
    ucsc_lc_df = length_checks_df(ucsc_lc_fpath,ucsc_taxid_dict)
    for ucsc_tid, row in any_ncbi_ids.iterrows():
        if not (recalc_pass_checks or ucsc_tid not in ucsc_lc_df.index):
            continue
        clade_tid_fpath = "{0}/{1}/{2}/{3}.fasta".format(reorg_parent,ds_id,clade_level,ucsc_tid)
        ucsc_df = fautil.load_UCSC_fasta_df(clade_tid_fpath)
        if length_how == "non_zero_clade":
            non_zero_ucsc = ucsc_df.loc[ucsc_df['length'] > 0, :]
            accepted_records = non_zero_ucsc.index
            accepted_med_len = np.median(ucsc_df.loc[accepted_records, 'length'])
        for taxid in ucsc_taxid_dict:
            check_col = "{0}_check".format(taxid)
            tax_tid_len = ucsc_df.loc[ucsc_df['NCBI_taxid']==taxid,'length'].iloc[0]
            ucsc_lc_df.loc[ucsc_tid,check_col] = tax_tid_len >= len_tol*accepted_med_len
    ucsc_lc_df.to_csv(ucsc_lc_fpath,sep='\t')
    return ucsc_lc_df

if __name__ == "__main__":
    pd.options.display.max_columns = None
    update_ncbi_checks = False
    if update_ncbi_checks:
        length_metrics_fpath = "{0}/length_metrics.tsv".format(dir_vars['combined_run'])
        # length_metric_checks(length_metrics_fpath,recalc_pass_checks=True,recalc_spec_checks=True)
        length_metric_checks(length_metrics_fpath, recalc_pass_checks=True, recalc_spec_checks=False)
    update_ucsc_checks = True
    if update_ucsc_checks:
        ucsc_lc_df = ucsc_length_checks()
    checks_fpath = "{0}/length_checks.tsv".format(dir_vars['combined_run'])
    check_df = pd.read_csv(checks_fpath, sep='\t', index_col=0)
    for col in check_df:
        col_nonna = check_df[col].dropna()
        print("Number of passed tests for label {0}: {1}".format(col,sum(col_nonna)))