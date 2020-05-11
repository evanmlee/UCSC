import os
import pandas as pd
import numpy as np
from IPython.display import display

from conservation import analysis_calc as ac, analysis_record_filter as ar_filt
from utility import UCSCerrors
from utility import fastaUtility as fasta



def gene_summary_table(align_df, test_idx,blos_df, unique_thresh=1,display_summary=False,
                       summary_table_outpath="", use_jsd_gap_penalty=True):
    """Given an align_df representing OrthoDB and NCBI record multiple sequence alignment, calculates JSD, BLOSUM,
    gap and variant metrics for the OrthoDB records only.

    :param align_df: OrthoDB and NCBI MSA DataFrame. Columns are 1-indexed alignment positions and index is record ids.
    :param test_idx: Record id of test_species (all other ODB records will be used as the outgroup for JSD and BLOSUM)
    :param blos_df: BLOSUM62 matrix as DataFrame taking amino acid chars as index/col values
    :param unique_thresh: Default 1 -> substitutions can be present in only test_species to be considered unique.
    If given as an int, will be interpreted as maximum number of instances for "uniqueness"; if given as a ratio,
    converted to corresponding number of instances scaled to number of sequences in align_df.
    :param display_summary: If true, displays calculated summary statistics table to stdout
    :param summary_table_outpath: If provided, will write summary table as tsv to path.
    :param: use_jsd_gap_penalty: If True, will use gap penalty in calculating
    :return:
    """
    #If unique thresh given as a ratio, convert to maximum number of instances to be considered unique. Else use as is.
    if type(unique_thresh) == float and unique_thresh < 1:
        unique_thresh = max(int(len(align_df) * unique_thresh), 1)
    uniques = ac.find_uniques(align_df,unique_thresh,test_idx)
    unique_pos = uniques.columns
    if len(unique_pos) == 0:
        raise UCSCerrors.SequenceAnalysisError(0,"No species unique substitutions under occurence " 
                                      "threshold {0} instances".format(unique_thresh))
    n_seq = len(align_df)
    #Calculate JSD and BLOSUM + z-scores for entire alignment
    jsd, jsd_z = ac.generate_jsd_series(test_idx,align_df,keep_test_spec=False,use_gap_penalty=use_jsd_gap_penalty)
    test_og_blos_srs, test_og_blos_z_srs = ac.test_outgroup_blosum_series(align_df,test_idx,blos_df)
    summary_col_labels = ['NCBI Species Position','NCBI Variant','NCBI Variant Count',
                          'Outgroup Variant', 'Outgroup Variant Count', 'Analysis Sequences', 'Gap Fraction',
                          'JSD','JSD Z-Score','Test-Outgroup BLOSUM62', 'Test-Outgroup BLOSUM Z-Score',
                          'Outgroup Pairwise BLOSUM62']
    summary_df = pd.DataFrame(columns=summary_col_labels)
    summary_df.index.name = "MSA Position"
    for pos in unique_pos:
        #Calculate variant count metrics, pull blos/jsd values from alignment-calculated Series
        aln_col = uniques.loc[:,pos]
        vc_metrics = ac.variant_counts(aln_col,test_idx)
        tv, tvc, ov, ovc, gf = [vc_metrics[label] for label in ['test_variant','test_variant_count','outgroup_variant',
                                                            'outgroup_variant_count','gap_fraction']]
        test_og_blos,test_og_blos_z = test_og_blos_srs[pos],test_og_blos_z_srs[pos]
        #Calculate outgroup pairwise BLOSUM, native sequence position of substitution in test-species
        og_pw_blos = ac.pairwise_outgroup_blosum(aln_col,test_idx,blos_df)
        native_pos = ac.align_pos_to_native_pos(align_df,test_idx,pos)
        row_dict = dict(zip(summary_col_labels,[native_pos,tv,tvc,ov,ovc,n_seq,gf,jsd[pos],
                                                jsd_z[pos],test_og_blos,test_og_blos_z,og_pw_blos]))
        summary_df.loc[pos,:] = row_dict
    if display_summary:
        print("Test Species Index: {0}".format(test_idx))
        display(summary_df)
    if summary_table_outpath:
        summary_df.to_csv(summary_table_outpath,sep='\t',float_format='%.5f')
    return summary_df

def load_summary_table(summary_fpath):
    summary_df = pd.read_csv(summary_fpath,sep='\t',index_col=0)
    return summary_df

def overall_summary_table(dir_vars, xref_table, taxid_dict,
                          tid_subset=[], UCSC_tax_subset=[],
                          use_jsd_gap_penalty=True,force_recalc=False,
                          modified_outpath="",length_checks_fpath="",skip_overall=True):
    """Calculates summary analysis statistics for every gene symbol in gene_symbols.

    :param (dict) dir_vars: Contains paths for run-specific directories, see directoryUtility.
    :param (DataFrame) xref_table: Maps UCSC transcript IDs to NCBI Gene IDs
    :param (dict) taxid_dict: maps NCBI taxids to species names
    :param (array-like) tid_subset: If provided, will only attempt analysis for transcript_ids in tid_subset

    :param (array-like) UCSC_tax_subset: If provided, limits UCSC sequence data to taxids in UCSC_tax_subset
    :param use_jsd_gap_penalty: Determines if JSD calculations are performed with gap penalty or not. If you change
    this, it is recommended to set force_recalc to True to avoid any inconsistencies between files calculated before
    the change and after.
    :param force_recalc: If True, recalculates and rewrites all summary statistic data.
    :param modified_outpath: If provided, writes overall summary table to modified_outpath.
    :param length_checks_fpath: If provided, reads booolean table from provided path and only runs gene_summary table
    for records which passed length checks (see analysis_record_filter.py)
    :return: None. Writes individual gene summary tables to appropriate output subdirectories and overall summary table
    to [run_name]/summary/overall_summary.tsv Overall_summary table contains new columns corresponding to 1) gene symbol,
    2) gene-specific MSA position, 3, 4) Unique Substitution Wide Z-scores for JSD and Test-Outgroup BLOSUM
    MSA position,
    """
    from utility.UCSCerrors import load_errors, write_errors, print_errors
    from utility import directoryUtility as dirutil
    from conservation.analysis_calc import blos_df
    summary_run_dir,bestNCBI_dir= dir_vars['summary_run'],dir_vars['bestNCBI_parent']
    errors_fpath = "{0}/errors.tsv".format(summary_run_dir)
    check_qes, query_error_df = load_errors(errors_fpath,error_type="NCBIQueryError")

    #Read length checks information if necessary
    if length_checks_fpath:
        length_checks = True
        lc_df = pd.read_csv(length_checks_fpath,sep='\t',index_col=0)
    else:
        length_checks = False
    #Columns for overall summary table
    overall_summary_col_labels = ['Transcript ID','Gene','MSA Position','NCBI Species Position', 'NCBI Variant',
                                  'NCBI Variant Count', 'Outgroup Variant', 'Outgroup Variant Count',
                                  'Analysis Sequences', 'Gap Fraction',
                                  'JSD', 'JSD Alignment Z-Score','JSD US Z-Score',
                                  'Test-Outgroup BLOSUM62', 'Test-Outgroup BLOSUM Alignment Z-Score',
                                  'Test-Outgroup BLOSUM US Z-Score','Outgroup Pairwise BLOSUM62']
    for ncbi_taxid in taxid_dict:
        taxid_parent_dir = "{0}/conservation/{1}".format(summary_run_dir,ncbi_taxid)
        overall_df = pd.DataFrame(columns=overall_summary_col_labels)
        if modified_outpath:
            overall_summary_fpath = "{0}/{1}_summary.tsv".format(modified_outpath, ncbi_taxid)
        else:
            overall_summary_fpath = "{0}/conservation/{1}_summary.tsv".format(summary_run_dir,ncbi_taxid)

        tax_analysis_errors_fpath = "{0}/{1}_analysis_errors.tsv".format(summary_run_dir,ncbi_taxid)
        check_aes, analysis_error_df = load_errors(tax_analysis_errors_fpath, error_type="SequenceAnalysisError")
        tax_record_errors_fpath = "{0}/{1}_record_errors.tsv".format(summary_run_dir, ncbi_taxid)
        check_res, record_error_df = load_errors(tax_record_errors_fpath)
        if length_checks:
            col_label = "{0}_check".format(ncbi_taxid)
            lc_col = lc_df[col_label]
        for tid,row in xref_table.iterrows():
            ncbi_hgid,hgnc_symbol = row['NCBI_gid'], row['HGNC_symbol']
            if modified_outpath:
                tid_subdir = "{0}/{1}/{2}".format(modified_outpath, ncbi_taxid,tid)
            else:
                tid_subdir = "{0}/{1}".format(taxid_parent_dir,tid)
            out_summary_fpath = "{0}/{1}_summary.tsv".format(tid_subdir, tid)
            out_records_fpath = "{0}/{1}_records.tsv".format(tid_subdir, tid)
            combined_fasta_fpath = "{0}/combined/{1}.fasta".format(bestNCBI_dir, tid)

            if not os.path.exists(out_summary_fpath) or force_recalc:
                if check_qes and tid in query_error_df['tid'].unique():
                    print("Transcript ID {0} failed to fetch NCBI orthologs".format(tid))
                    continue
                if check_aes and tid in analysis_error_df['tid'].unique():
                    print_errors(analysis_error_df,tid,error_type="SequenceAnalysisError")
                    continue
                #Check for tid_subset and length_checks
                if ((len(tid_subset) > 0 and tid in tid_subset) or len(tid_subset) == 0) \
                        and (length_checks and lc_col[tid]==True):
                    combined_fasta_df = fasta.load_UCSC_NCBI_df(combined_fasta_fpath,taxid_dict)
                    # Default behavior: Drop other NCBI records when calculating unique substitutions for individual
                    # NCBI species
                    ucsc_partititon = (combined_fasta_df['NCBI_taxid'].isin(UCSC_tax_subset)) & \
                                        (combined_fasta_df.index.str.contains("ENST"))
                    ncbi_partition = (combined_fasta_df['NCBI_taxid']==ncbi_taxid) & \
                                     ~(combined_fasta_df.index.str.contains("ENST"))
                    ucsc_analysis_subset = combined_fasta_df.loc[ucsc_partititon].index
                    ncbi_idx = combined_fasta_df.loc[ncbi_partition,:].index
                    #Filter analysis subset. Drop redundant UCSC data (ie UCSC data representing species in NCBI
                    #analaysis subset or close evolutionary relatives which are
                    try:
                        dirutil.create_directory(tid_subdir)
                        records_df, filtered_aln_path = ar_filt.filter_analysis_subset(combined_fasta_fpath,out_records_fpath,
                                                                UCSC_analysis_subset=ucsc_analysis_subset,NCBI_record_subset=ncbi_idx,
                                                                taxid_dict=taxid_dict,drop_redundant=True,
                                                                      drop_ncbi_from_ucsc=True)
                        align_df = fasta.align_fasta_to_df(filtered_aln_path,ucsc_flag=True)
                        summary_df = gene_summary_table(align_df,ncbi_idx,blos_df,
                                                    summary_table_outpath=out_summary_fpath,
                                                        use_jsd_gap_penalty=use_jsd_gap_penalty)
                    except UCSCerrors.SequenceAnalysisError as sae:
                        dirutil.remove_thing(tid_subdir)
                        write_errors(tax_analysis_errors_fpath,tid,sae)
                        continue
                else:
                    if check_res and tid in record_error_df['tid'].unique():
                        print_errors(record_error_df,tid)
                    elif length_checks and np.isnan(lc_col[tid]):
                        emsg = "No NCBI record data for taxid {0} present.".format(ncbi_taxid)
                        record_error = UCSCerrors.SequenceDataError(0,emsg)
                        write_errors(tax_record_errors_fpath,tid,record_error)
                    elif length_checks and lc_col[tid]==False:
                        emsg = "NCBI record data for taxid {0} failed record length quality checks.".format(ncbi_taxid)
                        record_error = UCSCerrors.SequenceDataError(1, emsg)
                        write_errors(tax_record_errors_fpath, tid, record_error)
                    continue
            else:
                summary_df = load_summary_table(out_summary_fpath)
            if not skip_overall:
                formatted = summary_df.reset_index(drop=False)
                formatted.insert(0,"Gene",[hgnc_symbol]*len(formatted))
                formatted.insert(0, "Transcript ID", [tid] * len(formatted))
                formatted = formatted.rename(columns={'JSD Z-Score':'JSD Alignment Z-Score','Test-Outgroup BLOSUM Z-Score':
                                                'Test-Outgroup BLOSUM Alignment Z-Score'})
                overall_df = overall_df.append(formatted, ignore_index=True, sort=False)
        # #Unique Substitution Set-wide Z scores for JSD and BLOSUM
        if not skip_overall:
            us_jsd, us_blos = ac.calc_z_scores(overall_df['JSD']), ac.calc_z_scores(overall_df['Test-Outgroup BLOSUM62'])
            overall_df.loc[:,'JSD US Z-Score'] = us_jsd
            overall_df.loc[:,'Test-Outgroup BLOSUM US Z-Score'] = us_blos
            overall_df.to_csv(overall_summary_fpath, sep='\t', float_format='%.5f')


